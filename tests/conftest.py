"""Small LAMMPS data files written on the fly, so no test depends on a shipped dataset."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


def write_lammps(
    path: Path,
    coordinates: np.ndarray,
    bonds: list[tuple[int, int]],
    box: tuple[float, float, float],
    lo: tuple[float, float, float] = (0.0, 0.0, 0.0),
    atom_ids: list[int] | None = None,
) -> Path:
    """A full-style data file; atom IDs default to 1..N and bonds are atom-ID pairs."""
    ids = atom_ids if atom_ids is not None else list(range(1, len(coordinates) + 1))
    with open(path, "w") as out:
        out.write("test network\n\n")
        out.write(f"{len(ids)} atoms\n1 atom types\n{len(bonds)} bonds\n1 bond types\n\n")
        for axis, name in enumerate("xyz"):
            out.write(f"{lo[axis]:.6f} {lo[axis] + box[axis]:.6f} {name}lo {name}hi\n")
        out.write("\nMasses\n\n1 1.0\n\nAtoms # full\n\n")
        for atom_id, (x, y, z) in zip(ids, coordinates):
            out.write(f"{atom_id} 1 1 0.000000 {x:.6f} {y:.6f} {z:.6f}\n")
        out.write("\nBonds\n\n")
        for bond_id, (a, b) in enumerate(bonds, start=1):
            out.write(f"{bond_id} 1 {a} {b}\n")
    return path


def random_periodic_network(n: int, box: tuple[float, float, float], k: int, seed: int):
    """Random points in a periodic box, each bonded to its k nearest neighbours (minimum image)."""
    rng = np.random.default_rng(seed)
    lengths = np.array(box)
    positions = rng.random((n, 3)) * lengths
    delta = positions[:, None, :] - positions[None, :, :]
    delta -= np.round(delta / lengths) * lengths
    distance = np.linalg.norm(delta, axis=-1)
    np.fill_diagonal(distance, np.inf)
    bonds = set()
    for i in range(n):
        for j in np.argsort(distance[i])[:k]:
            bonds.add((min(i, j) + 1, max(i, j) + 1))
    return positions, sorted(bonds)


@pytest.fixture
def ring_file(tmp_path: Path) -> Path:
    """12 atoms in a ring along x; the last bond crosses the periodic boundary."""
    coordinates = np.array([[0.5 + i, 5.0, 5.0] for i in range(12)])
    bonds = [(i, i + 1) for i in range(1, 12)] + [(12, 1)]
    return write_lammps(tmp_path / "ring.data", coordinates, bonds, box=(12.0, 10.0, 10.0))


@pytest.fixture
def chain_file(tmp_path: Path) -> Path:
    """An open chain that never reaches the periodic boundary: nothing wraps."""
    coordinates = np.array([[2.0 + i, 5.0, 5.0] for i in range(5)])
    bonds = [(i, i + 1) for i in range(1, 5)]
    return write_lammps(tmp_path / "chain.data", coordinates, bonds, box=(20.0, 10.0, 10.0))


@pytest.fixture
def random_file(tmp_path: Path) -> Path:
    """250 atoms, 468 bonds. Golden results in test_shortest_path.py were checked against the
    reference shortest-path code on exactly this network."""
    positions, bonds = random_periodic_network(250, (14.0, 10.0, 12.0), k=3, seed=1)
    return write_lammps(tmp_path / "random.data", positions, bonds, box=(14.0, 10.0, 12.0))
