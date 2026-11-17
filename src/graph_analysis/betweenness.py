"""Native spatially-weighted edge betweenness (GEBC/OBC/DBC) for LAMMPS polymer networks.

Computes geodesic edge betweenness centrality (GEBC, the standard topological
metric) alongside two spatially-weighted variants in a single pass: OBC
(orientation-weighted, by alignment with a loading axis) and DBC
(distance-weighted, by source-target separation).

Requires the in-house igraph build with `Graph.edge_betweenness_spatial` (see
`vendor/README.md` at the repo root for how to build it) - the stock PyPI
igraph package does not have this method, and importing this subpackage
without it raises a clear `RuntimeError` at analysis time, not at import time.

The input is a `lammps.LammpsNetwork` (atom IDs, coordinates, box bounds), not the
`networkx`-backed `network.PolymerNetwork` of the topology analyses, so this module is
deliberately not re-exported from the top-level `graph_analysis` package. Import it
explicitly: `from graph_analysis.betweenness import SpatialBetweennessAnalysis`.

One native, single-threaded GEBC/OBC/DBC traversal runs per analysis. OBC and DBC
weight each source-target pair by its minimum-image displacement, while the
shortest path between the two atoms on the periodic graph may end on another
periodic image of the target; this happens more often the longer the path is
compared with the box.
"""

from __future__ import annotations

import os
import pickle
import platform
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from numbers import Real
from pathlib import Path
from time import perf_counter

import igraph as ig
import numpy as np

from .lammps import Direction, FloatArray, LammpsDataReader, LammpsFormatError, LammpsNetwork

__all__ = [
    "AnalysisConfig",
    "CentralityResult",
    "CentralityValidationError",
    "Direction",
    "LammpsDataReader",
    "LammpsFormatError",
    "LammpsNetwork",
    "SpatialBetweennessAnalysis",
    "ValidationReport",
]


@dataclass(frozen=True, slots=True)
class AnalysisConfig:
    """Native analysis settings; dbc_length is used only for output postprocessing.

    None leaves DBC without any reference-length division. Raw DBC always uses
    source-target distance alone, regardless of this setting.
    """

    direction: Direction | str | int = Direction.X
    validate_stock: bool = False
    rtol: float = 1e-12
    atol: float = 1e-8
    dbc_length: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "direction", Direction.parse(self.direction))
        if self.dbc_length is not None:
            if isinstance(self.dbc_length, (bool, np.bool_)) or not isinstance(self.dbc_length, Real):
                raise ValueError("dbc_length must be a finite, positive number or None")
            try:
                length = float(self.dbc_length)
            except OverflowError as exc:
                raise ValueError("dbc_length must be finite and positive") from exc
            if not np.isfinite(length) or length <= 0:
                raise ValueError("dbc_length must be finite and positive")
            object.__setattr__(self, "dbc_length", length)
        for name in ("rtol", "atol"):
            value = getattr(self, name)
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")


@dataclass(frozen=True, slots=True)
class ValidationReport:
    max_absolute_difference: float
    mean_absolute_difference: float
    max_relative_difference: float
    pearson: float | None
    allclose: bool
    rtol: float
    atol: float
    relative_denominator_floor: float = 1e-15

    @classmethod
    def compare(cls, actual: FloatArray, reference: FloatArray, *, rtol: float, atol: float) -> ValidationReport:
        if actual.shape != reference.shape:
            raise ValueError("Centrality arrays must have identical shapes")
        if not np.isfinite(actual).all() or not np.isfinite(reference).all():
            raise ArithmeticError("Cannot validate nonfinite centrality values")
        difference = np.abs(actual - reference)
        pearson = None
        if actual.size > 1 and np.std(actual) > 0 and np.std(reference) > 0:
            pearson = float(np.corrcoef(actual, reference)[0, 1])
        return cls(float(difference.max(initial=0)), float(difference.mean()) if difference.size else 0.0,
                   float((difference / np.maximum(np.abs(reference), 1e-15)).max(initial=0)), pearson,
                   bool(np.allclose(actual, reference, rtol=rtol, atol=atol)), rtol, atol)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class CentralityResult:
    """Raw undirected GEBC/OBC/DBC arrays aligned with ``network.bond_ids``."""

    network: LammpsNetwork
    raw_gebc: FloatArray
    raw_obc: FloatArray
    raw_dbc: FloatArray
    neff: int
    metadata: dict[str, object]
    validation: ValidationReport | None = None

    def __post_init__(self) -> None:
        if not 0 <= self.neff <= self.network.n_atoms:
            raise ValueError("neff is outside the network's vertex range")
        for name in ("raw_gebc", "raw_obc", "raw_dbc"):
            values = np.array(getattr(self, name), dtype=np.float64, copy=True)
            if values.shape != (self.network.n_bonds,) or not np.isfinite(values).all() or np.any(values < 0):
                raise ArithmeticError(f"{name} must contain one finite, nonnegative value per bond")
            values.setflags(write=False)
            object.__setattr__(self, name, values)

    @property
    def normalization_factor(self) -> float:
        return 2.0 / (self.neff * (self.neff - 1)) if self.neff > 1 else 0.0

    @property
    def obc(self) -> FloatArray:
        """OBC normalized by the number of unordered pairs of nonisolated nodes."""
        return self.raw_obc * self.normalization_factor

    @property
    def dbc(self) -> FloatArray:
        """Pair-normalized DBC, optionally divided by the selected postprocessing length."""
        length = self.metadata.get("dbc_length")
        if length is not None:
            return self.dbc_for_length(length, normalized=True)
        return self.raw_dbc * self.normalization_factor

    def dbc_for_length(self, length: float, *, normalized: bool = False) -> FloatArray:
        """Apply any positive reference length to the same unmodified raw DBC array."""
        if length is None:
            raise ValueError("dbc_length must be finite and positive")
        length = AnalysisConfig(dbc_length=length).dbc_length
        with np.errstate(over="raise", invalid="raise"):
            values = self.raw_dbc / length
        return values * self.normalization_factor if normalized else values

    def summary(self) -> dict[str, object]:
        report = dict(self.metadata)
        report.update(neff=self.neff, normalization_factor=self.normalization_factor,
                      gebc_validation=self.validation.to_dict() if self.validation else None)
        return report

    def save(self, path: str | Path, *, overwrite: bool = False) -> Path:
        """Pickle this whole result to one file, via a temporary file before replacement.

        Matches the rest of this package (`ScanlanCaseResult`, `NetworkStats`, ...):
        one plain pickle of the dataclass, not a NumPy archive plus a side-car report.
        `result.network` carries `atom_ids`/`edges`/`bond_atom_ids`/box info; `.obc`/
        `.dbc` are computed properties, not stored separately. Existing output
        requires overwrite=True.
        """
        output = Path(path).expanduser().resolve()
        if output.suffix != ".pkl":
            raise ValueError("Output filename must end in .pkl")
        if self.network.source is not None and Path(self.network.source) == output:
            raise ValueError("Output file must not replace the input LAMMPS file")
        if not overwrite and output.exists():
            raise FileExistsError(f"Output already exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(dir=output.parent, prefix=f".{output.name}.", delete=False) as stream:
                temporary = Path(stream.name)
                pickle.dump(self, stream, protocol=pickle.HIGHEST_PROTOCOL)
            # NamedTemporaryFile creates the file private (0600); give it the mode open() would.
            umask = os.umask(0)
            os.umask(umask)
            os.chmod(temporary, 0o666 & ~umask)
            os.replace(temporary, output)
            temporary = None
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        return output


class CentralityValidationError(RuntimeError):
    """Custom GEBC disagrees with stock; spatial results must not be accepted."""

    def __init__(self, report: ValidationReport) -> None:
        self.report = report
        super().__init__(f"GEBC validation failed: {report.to_dict()}")


class SpatialBetweennessAnalysis:
    """Prepare a graph once and compute all three centrality channels together.

    The graph is undirected, simple and at bead/bond level. Shortest paths are
    always unweighted; geometry only weights source-target contributions.
    The library does not install global igraph progress handlers.
    """

    def __init__(self, network: LammpsNetwork, config: AnalysisConfig | None = None) -> None:
        self.network = network
        self.config = config or AnalysisConfig()
        self._graph = ig.Graph(n=network.n_atoms, edges=network.edges, directed=False)
        if not self._graph.is_simple():
            raise ValueError("Expected a simple atom/bond graph; found self-loops or duplicate atom-pair bonds")
        edge_ids = np.asarray(self._graph.get_edgelist(), dtype=np.int64).reshape(-1, 2)
        if not np.array_equal(edge_ids, np.sort(network.edges, axis=1)):
            raise RuntimeError("igraph edge IDs do not match LAMMPS bond-row order")
        self._neff = int(np.count_nonzero(self._graph.degree()))
        self._sample_parent: tuple[int, int] | None = None

    @classmethod
    def from_lammps(cls, path: str | Path, config: AnalysisConfig | None = None) -> SpatialBetweennessAnalysis:
        return cls(LammpsDataReader().read(path), config)

    @property
    def dbc_length(self) -> float | None:
        """Optional postprocessing denominator; never applied to raw DBC."""
        return self.config.dbc_length

    def sample(self, n_vertices: int) -> SpatialBetweennessAnalysis:
        """Return a BFS-neighborhood induced subgraph for debugging, not an estimate.

        The seed is the smallest atom ID. If its component is smaller than the
        requested size, all vertices reachable from that seed are used.
        """
        if isinstance(n_vertices, bool) or not isinstance(n_vertices, (int, np.integer)) or n_vertices < 1:
            raise ValueError("Sample size must be a positive integer")
        if self.network.n_atoms == 0:
            raise ValueError("Cannot sample an empty graph")
        vertices = np.asarray(self._graph.bfs(0)[0][:n_vertices], dtype=np.int64)
        result = SpatialBetweennessAnalysis(self.network.subset(vertices), self.config)
        result._sample_parent = (self.network.n_atoms, self.network.n_bonds)
        return result

    def inspect(self) -> dict[str, object]:
        """Input/topology metadata; no centrality calculation is performed."""
        return {
            "schema_version": 2,
            "raw_dbc_weight": "distance",
            "input": str(self.network.source) if self.network.source else None,
            "nodes": self.network.n_atoms,
            "edges": self.network.n_bonds,
            "neff": self._neff,
            "components": len(self._graph.connected_components()),
            "box_lengths": self.network.box_lengths.tolist(),
            "direction": int(self.config.direction),
            "loading_axis": self.config.direction.name.lower(),
            "dbc_length": self.dbc_length,
            "dbc_length_source": "none" if self.config.dbc_length is None else "explicit",
            "edge_order_verified": True,
            "graph_level": "atom/bond",
            "pair_convention": "unordered, raw undirected",
            "periodic_axes": ["x", "y", "z"],
            "sample": self._sample_parent is not None,
            "sample_parent_nodes_edges": self._sample_parent,
            "python": platform.python_version(),
            "platform": platform.platform(),
            "python_igraph": ig.__version__,
            "c_igraph": ig.__igraph_version__,
        }

    def run(self) -> CentralityResult:
        """Compute raw channels and optionally check GEBC before returning results."""
        native = getattr(self._graph, "edge_betweenness_spatial", None)
        if not callable(native):
            raise RuntimeError("The modified igraph build is required. Use vendor/.spatial-venv/bin/python, or build it first: cd vendor && bash python-igraph/scripts/build_spatial.sh (see vendor/README.md).")
        if getattr(ig._igraph, "SPATIAL_DBC_WEIGHT", None) != "distance":
            raise RuntimeError("The native DBC convention is outdated. Rebuild both libraries (cd vendor && bash python-igraph/scripts/build_spatial.sh) before using distance-only raw DBC.")
        metadata = self.inspect()
        metadata["started_at_utc"] = datetime.now(timezone.utc).isoformat()
        start = perf_counter()
        raw_gebc, raw_obc, raw_dbc = [np.asarray(values, dtype=np.float64) for values in native(self.network.coordinates, self.network.box_lengths, int(self.config.direction))]
        metadata["spatial_seconds"] = perf_counter() - start
        validation = None
        if self.config.validate_stock:
            start = perf_counter()
            stock = np.asarray(self._graph.edge_betweenness(directed=False), dtype=np.float64)
            metadata["stock_seconds"] = perf_counter() - start
            validation = ValidationReport.compare(raw_gebc, stock, rtol=self.config.rtol, atol=self.config.atol)
            if not validation.allclose:
                raise CentralityValidationError(validation)
        result = CentralityResult(self.network, raw_gebc, raw_obc, raw_dbc, self._neff, metadata, validation)
        try:
            import resource
        except ImportError:
            metadata["peak_process_rss_mib"] = None
        else:
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            metadata["peak_process_rss_mib"] = rss / (1024**2 if platform.system() == "Darwin" else 1024)
        return result
