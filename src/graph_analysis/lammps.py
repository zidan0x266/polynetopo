"""LAMMPS data input shared by the spatial analyses (betweenness, shortest paths).

`LammpsDataReader` streams a data file (atom style full, molecular, bond or angle) into a
`LammpsNetwork` (atom IDs, coordinates, bonds, orthorhombic box) and validates it once, in one
place. `Direction`
is the axis that betweenness loads along and shortest paths wrap along.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from operator import index
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


class Direction(IntEnum):
    """Loading axis; values match the native igraph API."""

    X = 0
    Y = 1
    Z = 2

    @classmethod
    def parse(cls, value: Direction | str | int) -> Direction:
        if isinstance(value, str):
            value = {"x": 0, "y": 1, "z": 2, "0": 0, "1": 1, "2": 2}.get(value.strip().lower(), -1)
        if isinstance(value, (bool, np.bool_)):
            raise ValueError("Direction must be x, y, z, 0, 1 or 2")
        try:
            return cls(index(value))
        except (TypeError, ValueError) as exc:
            raise ValueError("Direction must be x, y, z, 0, 1 or 2") from exc


def _integers(values: IntArray, name: str) -> IntArray:
    array = np.asarray(values)
    if array.size and not np.issubdtype(array.dtype, np.integer):
        raise ValueError(f"{name} must contain integers")
    if array.size and np.issubdtype(array.dtype, np.unsignedinteger) and array.max() > np.iinfo(np.int64).max:
        raise ValueError(f"{name} exceeds the int64 range")
    return np.array(array, dtype=np.int64, copy=True)


@dataclass(frozen=True, slots=True)
class LammpsNetwork:
    """An orthorhombic bead-level graph read from a LAMMPS data file, with explicit ID mapping.

    Distinct from `network.PolymerNetwork` (a `networkx` graph with no coordinates,
    used by the topology analyses): this one carries coordinates and the box, which
    the spatial analyses (betweenness, shortest paths) need.


    Coordinate rows follow sorted ``atom_ids``. ``edges`` contains zero-based
    indices into those rows. Bond arrays retain their original file order.
    Arrays are copied and made read-only so analysis cannot silently alter input.
    ``source`` is the resolved input path as a string, so that pickled results
    load under any Python version.
    """

    atom_ids: IntArray
    coordinates: FloatArray
    bond_ids: IntArray
    bond_types: IntArray
    edges: IntArray
    box_bounds: FloatArray
    source: str | None = None

    def __post_init__(self) -> None:
        for name in ("atom_ids", "bond_ids", "bond_types", "edges"):
            object.__setattr__(self, name, _integers(getattr(self, name), name))
        for name in ("coordinates", "box_bounds"):
            object.__setattr__(self, name, np.array(getattr(self, name), dtype=np.float64, copy=True))
        if self.atom_ids.ndim != 1 or self.bond_ids.ndim != 1:
            raise ValueError("Atom IDs and bond IDs must be one-dimensional")
        if self.coordinates.shape != (self.n_atoms, 3) or not np.isfinite(self.coordinates).all():
            raise ValueError("Coordinates must be a finite (n_atoms, 3) array")
        if np.any(self.atom_ids <= 0) or np.any(np.diff(self.atom_ids) <= 0):
            raise ValueError("Atom IDs must be positive, unique and sorted")
        if np.any(self.bond_ids <= 0) or np.unique(self.bond_ids).size != self.n_bonds:
            raise ValueError("Bond IDs must be positive and unique")
        if self.bond_types.shape != (self.n_bonds,) or np.any(self.bond_types <= 0):
            raise ValueError("Each bond must have a positive bond type")
        if self.edges.shape != (self.n_bonds, 2):
            raise ValueError("Edges must have shape (n_bonds, 2)")
        if self.edges.size and (self.edges.min() < 0 or self.edges.max() >= self.n_atoms):
            raise ValueError("Edge endpoint outside the vertex-index range")
        if self.box_bounds.shape != (3, 2) or not np.isfinite(self.box_bounds).all():
            raise ValueError("Box bounds must be a finite (3, 2) array")
        if not np.isfinite(self.box_lengths).all() or np.any(self.box_lengths <= 0):
            raise ValueError("Box lengths must be finite and strictly positive")
        if self.source is not None:
            object.__setattr__(self, "source", str(Path(self.source).resolve()))
        for name in ("atom_ids", "coordinates", "bond_ids", "bond_types", "edges", "box_bounds"):
            getattr(self, name).setflags(write=False)

    @property
    def n_atoms(self) -> int:
        return self.atom_ids.size

    @property
    def n_bonds(self) -> int:
        return self.bond_ids.size

    @property
    def box_lengths(self) -> FloatArray:
        return self.box_bounds[:, 1] - self.box_bounds[:, 0]

    @property
    def bond_atom_ids(self) -> IntArray:
        """Original LAMMPS endpoint IDs, with endpoint and bond-row order preserved."""
        return self.atom_ids[self.edges]

    def subset(self, vertices: IntArray) -> LammpsNetwork:
        """Build an induced debug subgraph; its centralities differ from the full graph."""
        vertices = np.unique(_integers(vertices, "vertices"))
        if vertices.size and (vertices.min() < 0 or vertices.max() >= self.n_atoms):
            raise ValueError("Sample vertex index outside the graph")
        selected = np.zeros(self.n_atoms, dtype=bool)
        selected[vertices] = True
        keep = selected[self.edges].all(axis=1)
        edges = np.searchsorted(vertices, self.edges[keep])
        return LammpsNetwork(self.atom_ids[vertices], self.coordinates[vertices], self.bond_ids[keep],
                              self.bond_types[keep], edges, self.box_bounds, self.source)


class LammpsFormatError(ValueError):
    """Invalid or unsupported LAMMPS data, with file/line context where possible."""


# Atom styles the reader accepts: the column of x in an Atoms row (id mol type [q] x y z [ix iy iz]).
ATOM_STYLES = {"full": 4, "molecular": 3, "bond": 3, "angle": 3}


class LammpsDataReader:
    """Read ``Atoms`` and ``Bonds`` from a 3D orthorhombic data file.

    The Atoms section may be of style full (``id mol type q x y z``), molecular, bond or angle
    (``id mol type x y z``), each with optional image flags; a bare ``Atoms`` header is
    interpreted as full style. Image flags are ignored
    because the native calculation uses minimum-image source-target distances.
    Atom IDs may be nonconsecutive. Unrelated numeric sections are skipped.
    """

    def read(self, path: str | Path) -> LammpsNetwork:
        source = Path(path).expanduser().resolve()
        counts: dict[str, int] = {}
        bounds: dict[str, tuple[float, float]] = {}
        atom_ids = coordinates = bonds = None
        x_column = ATOM_STYLES["full"]
        nread_atoms = nread_bonds = 0
        section = None
        in_header = True
        with source.open(encoding="utf-8") as stream:
            next(stream, None)  # LAMMPS reserves the first line for the title.
            for line_number, line in enumerate(stream, start=2):
                data, _, comment = line.partition("#")
                fields = data.split()
                if not fields:
                    continue
                try:
                    if fields[-3:] == ["xy", "xz", "yz"] or fields[-1] in ("avec", "bvec", "cvec"):
                        raise ValueError("Triclinic cells are unsupported; an orthorhombic box is required")
                    if in_header and len(fields) == 2 and fields[1] in ("atoms", "bonds"):
                        key, count = fields[1], int(fields[0])
                        if count < 0 or key in counts:
                            raise ValueError(f"Invalid or repeated {key} count")
                        counts[key] = count
                        continue
                    if in_header and len(fields) == 4 and fields[2:] in (["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]):
                        axis = fields[2][0]
                        if axis in bounds:
                            raise ValueError(f"Repeated {axis} box bounds")
                        bounds[axis] = (float(fields[0]), float(fields[1]))
                        continue
                    if fields[0][0].isalpha():
                        in_header = False
                        section = fields[0]
                        if section == "Atoms":
                            if "atoms" not in counts or atom_ids is not None:
                                raise ValueError("Missing atom count or repeated Atoms section")
                            style = comment.split()[0] if comment.split() else "full"
                            if style not in ATOM_STYLES:
                                raise ValueError(f"Unsupported atom style {style!r}; expected one of "
                                                 f"{', '.join(ATOM_STYLES)}")
                            x_column = ATOM_STYLES[style]
                            atom_ids = np.empty(counts["atoms"], dtype=np.int64)
                            coordinates = np.empty((counts["atoms"], 3), dtype=np.float64)
                        elif section == "Bonds":
                            if bonds is not None:
                                raise ValueError("Repeated Bonds section")
                            bonds = np.empty((counts.get("bonds", 0), 4), dtype=np.int64)
                        continue
                    if section == "Atoms":
                        if len(fields) not in (x_column + 3, x_column + 6) or nread_atoms >= len(atom_ids):
                            raise ValueError(f"Expected exactly the declared number of atom rows ({x_column + 3} or "
                                             f"{x_column + 6} columns for this atom style)")
                        atom_ids[nread_atoms] = int(fields[0])
                        coordinates[nread_atoms] = [float(value) for value in fields[x_column:x_column + 3]]
                        nread_atoms += 1
                    elif section == "Bonds":
                        if len(fields) != 4 or nread_bonds >= len(bonds):
                            raise ValueError("Expected exactly the declared number of bond rows (4 columns)")
                        bonds[nread_bonds] = [int(value) for value in fields]
                        nread_bonds += 1
                except (ValueError, OverflowError) as exc:
                    raise LammpsFormatError(f"{source}:{line_number}: {exc}") from exc

        if "atoms" not in counts or nread_atoms != counts["atoms"] or nread_bonds != counts.get("bonds", 0):
            raise LammpsFormatError(f"{source}: Atoms/Bonds sections do not match the declared counts")
        if set(bounds) != {"x", "y", "z"}:
            raise LammpsFormatError(f"{source}: missing orthorhombic box bounds")
        if atom_ids is None:
            atom_ids = np.empty(0, dtype=np.int64)
            coordinates = np.empty((0, 3), dtype=np.float64)
        if bonds is None:
            bonds = np.empty((0, 4), dtype=np.int64)
        order = np.argsort(atom_ids)
        atom_ids, coordinates = atom_ids[order], coordinates[order]
        if np.any(atom_ids <= 0) or np.any(np.diff(atom_ids) <= 0):
            raise LammpsFormatError(f"{source}: atom IDs must be positive and unique")
        edges = np.searchsorted(atom_ids, bonds[:, 2:4])
        if edges.size and (np.any(edges >= atom_ids.size) or not np.array_equal(atom_ids[edges], bonds[:, 2:4])):
            raise LammpsFormatError(f"{source}: a bond references an unknown atom ID")
        try:
            return LammpsNetwork(atom_ids, coordinates, bonds[:, 0], bonds[:, 1], edges,
                                  np.array([bounds[axis] for axis in "xyz"]), source)
        except ValueError as exc:
            raise LammpsFormatError(f"{source}: {exc}") from exc
