"""Shortest paths across periodic boundaries, which govern bond rupture.

Implements the method of Yu & Jackson, "Shortest Paths Govern Bond Rupture in
Thermoset Networks" (Macromolecules, 2025), for a network under periodic
boundary conditions (PBC). A shortest path (SP) that wraps through the box can
be longer than one image, so the search runs on a series of multi-image graphs.
For each period p = 1, 2, ..., M:

1. Build G_p from p + 1 periodic images of the network (image 0 is the box as
   given), plus `buffer_images` extra images on each side. Bonds are kept
   inside each image, and bonds that cross the periodic boundary connect
   neighbouring images instead (never the outer boundary of the outermost
   image). The buffer images let a path leave the stack between its two ends,
   instead of being forced to stay inside it.
2. For every atom i, find the shortest path from atom i in image 0 to the same
   atom in image p, and normalize its length by p: D^g = (path length) / p, the
   edges per box length. Normalizing makes the results of different p
   comparable. A path is valid only if the start/end atom is the only atom it
   visits twice; otherwise it is discarded.
3. Delete redundant paths. Case 1: the same path found from different start
   atoms. Case 2: a path at larger p that repeats one at smaller p. Both are
   removed when an already kept path of equal D^g contains all of a path's
   atoms (and Case 2 mostly by the validity rule above). Case 3 is optional
   (`min_sp_only`): among the SPs of one atom, keep only the shortest.

The search stops at the first p with no valid path, as the reference does: in
a cross-linked network, the paths at larger p then only repeat shorter ones.
This is a stopping rule rather than a guarantee; a structure that winds through
the box more than once before it closes, such as a single closed chain that
wraps the box twice, has its first valid path at a larger p. The first
(shortest) SP predicts where the first bond breaks and, through D^g_min / L0,
the strain at which it breaks.

Port notes - behaviour matches the reference implementation (the Jackson lab's
shortest_path_thermoset code as run on the cluster, with `buffer_images`)
except where listed:

- Input comes from `lammps.LammpsDataReader`, so the box, coordinates and
  bonds are validated once, in one place, and atom IDs may be nonconsecutive.
  Paths are reported in *atom IDs* (raw LAMMPS IDs, like every other output of
  this package), not zero-based indices, and images are counted from image 0
  (the reference counts from the first buffer image).
- A period is given as p, as in the paper (p = 1..M). The reference takes the
  number of images n = p + 1 and tries `range(image_min, image_max)`.
- A bond crosses the periodic boundary when its endpoints are more than half a
  box length apart along the PBC axis, tested, as in the reference, on the
  coordinates as written; LAMMPS `write_data` writes them wrapped into the box.
- Each multi-image graph is built once, vectorized, in the reference's exact
  edge order (tie-breaking between equally short paths depends on it), instead
  of once per source atom in every worker.
- A worker keeps only valid paths, so the parent never holds one path per atom.
- Redundancy removal groups kept paths by D^g instead of scanning every kept
  path for every new one (quadratic in the number of paths).
- Parallel bonds (the same atom pair bonded twice) are an error: a path lists
  atoms, so it cannot say which of two parallel bonds it used.

Cost: one breadth-first search per source atom over a graph with p + 1 + 2 *
`buffer_images` copies of the network, so cost grows as atoms x images x (atoms
+ bonds), and memory with the images and with p (a path has about D^g * p
atoms). For a large all-atom network, reduce it to crosslinkers first (the
paper does this), or limit the sources with `slice_max`.
"""

from __future__ import annotations

import hashlib
import logging
import multiprocessing as mp
import os
import pickle
import tempfile
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import igraph as ig
import numpy as np
from tqdm import tqdm

from .lammps import Direction, LammpsDataReader, LammpsNetwork

logger = logging.getLogger(__name__)

IntArray = np.ndarray
FloatArray = np.ndarray


@dataclass(frozen=True, slots=True)
class ShortestPathConfig:
    """Settings for a shortest-path search.

    `max_period` is the paper's M: periods p = 1..M are searched, each wrapping
    p box lengths. `buffer_images` adds that many extra images on each side of
    the p + 1 so a path may wander outside them; 0 keeps every path inside.
    `slice_max` restricts the source atoms to those within that distance of the
    box's low face along `direction` (None: every atom); it reduces cost and
    changes the set of paths, since only the SPs found from these atoms are
    kept, so it is meant for quick trials.
    `remove_redundant` applies Cases 1 and 2; `min_sp_only` also applies the
    paper's optional Case 3 (only the shortest SP of each atom).
    """

    direction: Direction | str | int = Direction.X
    max_period: int = 8
    buffer_images: int = 1
    slice_max: Optional[float] = None
    remove_redundant: bool = True
    min_sp_only: bool = False
    ncpus: Optional[int] = None
    checkpoint_file: Optional[Path] = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "direction", Direction.parse(self.direction))
        if self.max_period < 1:
            raise ValueError("max_period must be at least 1")
        if self.buffer_images < 0:
            raise ValueError("buffer_images must not be negative")
        if self.slice_max is not None and self.slice_max <= 0:
            raise ValueError("slice_max must be positive")
        if self.ncpus is not None and self.ncpus < 1:
            raise ValueError("ncpus must be at least 1")
        if self.checkpoint_file is not None:
            object.__setattr__(self, "checkpoint_file", Path(self.checkpoint_file))


@dataclass(frozen=True, slots=True)
class ShortestPathResult:
    """The shortest paths that wrap through the periodic box along one axis.

    `lengths[k]` is the normalized length D^g of `paths[k]`, an (n, 2) integer
    array of (atom ID, image) rows from the atom in image 0 to itself in image
    p; an image of -1, or above p, is a buffer image the path passed through.
    `min_length_by_period` has one (p, minimum D^g) row for each p that had a
    valid path, found before any redundant path was removed.
    """

    lengths: FloatArray
    paths: list[IntArray]
    min_length_by_period: FloatArray
    direction: Direction
    box_length: float
    max_period: int = 8
    buffer_images: int = 1
    min_sp_only: bool = False
    source: Optional[str] = None

    @property
    def min_length(self) -> float:
        """D^g_min, the length of the first (shortest) SP."""
        return float(self.lengths.min())

    def summary(self) -> dict[str, object]:
        return {
            "paths": len(self.paths),
            "direction": self.direction.name.lower(),
            "box_length": self.box_length,
            "max_period": self.max_period,
            "buffer_images": self.buffer_images,
            "min_sp_only": self.min_sp_only,
            "min_length": self.min_length,
            "max_length": float(self.lengths.max()),
            "mean_length": float(self.lengths.mean()),
            "median_length": float(np.median(self.lengths)),
            "periods_with_paths": [int(p) for p in self.min_length_by_period[:, 0]],
        }

    def save(self, path: str | Path, *, overwrite: bool = False) -> Path:
        """Pickle this whole result to one file, via a temporary file before replacement."""
        output = Path(path).expanduser().resolve()
        if output.suffix != ".pkl":
            raise ValueError("Output filename must end in .pkl")
        if self.source is not None and Path(self.source) == output:
            raise ValueError("Output file must not replace the input LAMMPS file")
        if not overwrite and output.exists():
            raise FileExistsError(f"Output already exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        _atomic_pickle(self, output)
        return output


class PeriodicShortestPathAnalysis:
    """Finds the shortest paths across the periodic boundary of one network."""

    def __init__(self, network: LammpsNetwork, config: ShortestPathConfig | None = None) -> None:
        _require_simple(network)
        self.network = network
        self.config = config or ShortestPathConfig()

    @classmethod
    def from_lammps(
        cls, path: str | Path, config: ShortestPathConfig | None = None
    ) -> PeriodicShortestPathAnalysis:
        return cls(LammpsDataReader().read(path), config)

    def run(self) -> ShortestPathResult:
        config = self.config
        axis = int(config.direction)
        box_length = float(self.network.box_lengths[axis])
        # Relative to the box's low face, as every box length below is.
        along = self.network.coordinates[:, axis] - self.network.box_bounds[axis, 0]
        slice_max = box_length if config.slice_max is None else config.slice_max
        sources = np.flatnonzero(along < slice_max)

        params = self._checkpoint_params(slice_max)
        state = self._load_checkpoint(params)
        lengths, paths, by_period = (
            ([], [], []) if state is None else (state["lengths"], state["paths"], state["by_period"])
        )
        first_period = 1 if state is None else state["period"] + 1
        if state is not None:
            logger.info("Resuming from a checkpoint after period %d", state["period"])

        for period in range(first_period, config.max_period + 1):
            n_images = period + 1
            graph = _image_graph(
                self.network.edges, along, self.network.n_atoms, box_length, n_images + 2 * config.buffer_images
            )
            valid = self._search(graph, sources, period)
            logger.info("period %d: %d valid paths from %d sources", period, len(valid), sources.size)
            if not valid:
                break
            normalized = [(len(path) - 1) / period for path in valid]
            lengths.extend(normalized)
            paths.extend(valid)
            by_period.append((period, min(normalized)))
            self._save_checkpoint(params, period, lengths, paths, by_period)

        if not lengths:
            raise ValueError("No shortest path found; is the network connected across the boundary?")
        if config.remove_redundant:
            lengths, paths = _remove_redundant(lengths, paths)
        if config.min_sp_only:
            lengths, paths = _keep_shortest_per_atom(lengths, paths)
        logger.info("Found %d shortest paths", len(paths))

        if config.checkpoint_file is not None:
            config.checkpoint_file.unlink(missing_ok=True)
        return ShortestPathResult(
            lengths=np.array(lengths, dtype=np.float64),
            paths=paths,
            min_length_by_period=np.array(by_period, dtype=np.float64),
            direction=config.direction,
            box_length=box_length,
            max_period=config.max_period,
            buffer_images=config.buffer_images,
            min_sp_only=config.min_sp_only,
            source=self.network.source,
        )

    def _search(self, graph: ig.Graph, sources: IntArray, period: int) -> list[IntArray]:
        """Valid shortest paths from each source in image 0 to itself in image `period`."""
        n_atoms, buffer = self.network.n_atoms, self.config.buffer_images
        first = buffer * n_atoms  # vertex index offset of image 0 (after the buffer images)
        tasks = [(first + int(source), first + int(source) + period * n_atoms) for source in sources]
        label = f"Finding shortest paths (period {period})"

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # igraph warns for every unreachable target
            if self.config.ncpus is None or self.config.ncpus == 1 or len(tasks) < 2:
                _init_worker(graph, self.network.atom_ids, buffer)
                results = [_trace(task) for task in tqdm(tasks, desc=label, leave=False)]
            else:
                with mp.Pool(self.config.ncpus, _init_worker, (graph, self.network.atom_ids, buffer)) as pool:
                    results = list(tqdm(pool.imap(_trace, tasks, chunksize=64), total=len(tasks), desc=label, leave=False))
        return [path for path in results if path is not None]

    # --- checkpointing: lets a long search resume after a crash, one period at a time ---

    def _checkpoint_params(self, slice_max: float) -> dict[str, object]:
        """What makes saved paths valid for this run: the network itself (a digest of its atoms,
        coordinates, bonds and box, so an edited input file is not mistaken for the old one) and the
        search settings. max_period and the redundancy options are not included, so an interrupted run
        can be resumed with a different max_period or Cases; the checkpoint is deleted once the search
        completes."""
        digest = hashlib.sha256()
        for array in (self.network.atom_ids, self.network.coordinates, self.network.edges, self.network.box_bounds):
            digest.update(np.ascontiguousarray(array).tobytes())
        return {
            "network": digest.hexdigest(),
            "direction": int(self.config.direction),
            "buffer_images": self.config.buffer_images,
            "slice_max": slice_max,
        }

    def _load_checkpoint(self, params: dict[str, object]) -> Optional[dict[str, object]]:
        file = self.config.checkpoint_file
        if file is None or not file.exists():
            return None
        with open(file, "rb") as inp:
            state = pickle.load(inp)
        if state["params"] != params:
            logger.warning("Ignoring checkpoint %s: its settings differ from this run", file)
            return None
        return state

    def _save_checkpoint(self, params, period, lengths, paths, by_period) -> None:
        file = self.config.checkpoint_file
        if file is not None:
            state = {"params": params, "period": period, "lengths": lengths, "paths": paths, "by_period": by_period}
            _atomic_pickle(state, file)


# --- input check ---------------------------------------------------------------------


def _require_simple(network: LammpsNetwork) -> None:
    """Reject a repeated bond or a bond from an atom to itself."""
    edges = network.edges
    if edges.size == 0:
        return
    loops = edges[:, 0] == edges[:, 1]
    if loops.any():
        raise ValueError(f"Atom {int(network.atom_ids[edges[loops][0, 0]])} is bonded to itself")
    pairs, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    repeated = counts > 1
    if repeated.any():
        a, b = network.atom_ids[pairs[repeated][0]]
        raise ValueError(
            f"{int(repeated.sum())} atom pairs are bonded more than once (first: {int(a)} and {int(b)}); "
            "a path lists atoms, so it cannot tell parallel bonds apart"
        )


# --- graph construction and search -------------------------------------------------


def _image_graph(edges: IntArray, along: FloatArray, n_atoms: int, box_length: float, n_images: int) -> ig.Graph:
    """The multi-image graph of `n_images` images: vertex (atom a, image i) is i * n_atoms + a.

    Edges are added in the reference implementation's order - every in-image
    bond (image by image), then for each pair of neighbouring images each
    boundary-crossing bond's forward edge followed by its reverse edge - because
    igraph breaks ties between equally short paths by edge order.
    """
    u, v = edges[:, 0], edges[:, 1]
    half = box_length / 2
    inside, between = [], []

    for image in range(n_images):
        shift = image * box_length
        keep = np.abs((along[u] + shift) - (along[v] + shift)) <= half
        inside.append(np.column_stack((u[keep], v[keep])) + image * n_atoms)

    for image in range(n_images - 1):
        # (u in this image, v in the next) and (v in this image, u in the next): the
        # bond crosses the boundary in one of the two directions, never both.
        forward = (along[v] + (image + 1) * box_length) - (along[u] + image * box_length) < half
        reverse = (along[u] + (image + 1) * box_length) - (along[v] + image * box_length) < half
        candidates = np.empty((len(u), 2, 2), dtype=np.int64)
        candidates[:, 0, 0] = u + image * n_atoms
        candidates[:, 0, 1] = v + (image + 1) * n_atoms
        candidates[:, 1, 0] = v + image * n_atoms
        candidates[:, 1, 1] = u + (image + 1) * n_atoms
        between.append(candidates[np.stack((forward, reverse), axis=1)])

    all_edges = np.concatenate(inside + between) if len(u) else np.empty((0, 2), dtype=np.int64)
    return ig.Graph(n=n_images * n_atoms, edges=all_edges.tolist(), directed=False)


_worker_graph: Optional[ig.Graph] = None
_worker_atom_ids: Optional[IntArray] = None
_worker_buffer: int = 0


def _init_worker(graph: ig.Graph, atom_ids: IntArray, buffer: int) -> None:
    global _worker_graph, _worker_atom_ids, _worker_buffer
    _worker_graph, _worker_atom_ids, _worker_buffer = graph, atom_ids, buffer


def _trace(task: tuple[int, int]) -> Optional[IntArray]:
    """A valid shortest path as (atom ID, image) rows, or None if it is unreachable or invalid.

    Images are counted from image 0; vertex v is atom v % n_atoms in image v // n_atoms.
    """
    source, target = task
    vertices = _worker_graph.get_shortest_paths(source, to=target, output="vpath")[0]
    if len(vertices) < 2:
        return None
    ids = np.asarray(vertices, dtype=np.int64)
    n_atoms = _worker_atom_ids.size
    atoms = ids % n_atoms
    if not _is_valid(atoms):
        return None
    return np.column_stack((_worker_atom_ids[atoms], ids // n_atoms - _worker_buffer))


def _is_valid(atoms: IntArray) -> bool:
    """True when the only atom visited twice is the start/end atom."""
    return atoms.size - np.unique(atoms).size == 1


def _remove_redundant(lengths: list[float], paths: list[IntArray]) -> tuple[list[float], list[IntArray]]:
    """Drop a path when an already kept path of equal D^g contains all of its atoms.

    That removes identical paths found from different start atoms (Case 1) and
    paths at larger p that only repeat one at smaller p (Case 2). Equal D^g is
    exact float equality, as in the reference.
    """
    kept_lengths: list[float] = []
    kept_paths: list[IntArray] = []
    atoms_by_length: dict[float, list[frozenset[int]]] = {}
    for length, path in zip(lengths, paths):
        atoms = frozenset(path[:-1, 0].tolist())
        kept = atoms_by_length.setdefault(length, [])
        if any(atoms <= other for other in kept):
            continue
        kept.append(atoms)
        kept_lengths.append(length)
        kept_paths.append(path)
    return kept_lengths, kept_paths


def _keep_shortest_per_atom(lengths: list[float], paths: list[IntArray]) -> tuple[list[float], list[IntArray]]:
    """Case 3 (optional in the paper): of the paths starting at one atom keep the shortest.

    The first one wins a tie, and the result is ordered by start atom, as in the
    reference's commented-out `get_SP_new`.
    """
    best: dict[int, int] = {}
    for index, (length, path) in enumerate(zip(lengths, paths)):
        atom = int(path[0, 0])
        if atom not in best or length < lengths[best[atom]]:
            best[atom] = index
    order = [best[atom] for atom in sorted(best)]
    return [lengths[i] for i in order], [paths[i] for i in order]


def _atomic_pickle(obj: object, output: Path) -> None:
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(dir=output.parent, prefix=f".{output.name}.", delete=False) as stream:
            temporary = Path(stream.name)
            pickle.dump(obj, stream, protocol=pickle.HIGHEST_PROTOCOL)
        # NamedTemporaryFile creates the file private (0600); give it the mode open() would.
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(temporary, 0o666 & ~umask)
        os.replace(temporary, output)
    except BaseException:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise
