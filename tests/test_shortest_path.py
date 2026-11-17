"""Shortest paths across the periodic boundary, on a cubic lattice whose shortest paths are known."""

from __future__ import annotations

import os
import pickle

import numpy as np
import pytest

from graph_analysis.lammps import LammpsNetwork
from graph_analysis.shortest_path import PeriodicShortestPathAnalysis, ShortestPathConfig

N = 4  # lattice sites per box length


def lattice(shift: float = 0.0, extra_bonds=()) -> LammpsNetwork:
    """N^3 beads, one per unit cell of a periodic box of side N, each bonded to its six neighbours.

    Along x, the shortest path that wraps the box is a straight line of N bonds, one for each
    of the N^2 lines of beads; at p = 2 the only paths of 2N bonds run twice along one line and
    are invalid, so the search stops there.
    """
    site = lambda x, y, z: (x % N) * N * N + (y % N) * N + (z % N)
    coordinates = np.array([(x + 0.5, y + 0.5, z + 0.5) for x in range(N) for y in range(N) for z in range(N)])
    coordinates[0, 1] += shift
    edges = [(site(x, y, z), site(x + dx, y + dy, z + dz))
             for x in range(N) for y in range(N) for z in range(N) for dx, dy, dz in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    edges += list(extra_bonds)
    return LammpsNetwork(np.arange(1, N**3 + 1), coordinates, np.arange(1, len(edges) + 1), np.ones(len(edges), int),
                         np.array(edges), np.array([[0.0, N]] * 3))


def test_lattice_paths():
    result = PeriodicShortestPathAnalysis(lattice(), ShortestPathConfig(max_period=3)).run()
    assert result.min_length == N
    assert np.all(result.lengths == N)
    assert len(result.paths) == N * N  # one per line of beads along x, after removing the repeats (Case 1)
    for path in result.paths:
        rows = path[:, 0] - 1
        assert len({row % (N * N) for row in rows}) == 1  # every bead of the path has the same y and z
        assert path[0, 0] == path[-1, 0] and (path[0, 1], path[-1, 1]) == (0, 1)
    assert result.min_length_by_period.tolist() == [[1.0, N]]


def test_min_sp_only_keeps_one_path_per_start_atom():
    result = PeriodicShortestPathAnalysis(lattice(), ShortestPathConfig(max_period=3, min_sp_only=True)).run()
    starts = [int(path[0, 0]) for path in result.paths]
    assert starts == sorted(set(starts))


def test_repeated_bond_is_rejected():
    with pytest.raises(ValueError, match="bonded more than once"):
        PeriodicShortestPathAnalysis(lattice(extra_bonds=[(0, 1)]))


def test_checkpoint_resumes_only_the_same_network(tmp_path):
    checkpoint = tmp_path / "checkpoint.pkl"
    config = ShortestPathConfig(max_period=3, checkpoint_file=checkpoint)
    interrupted = PeriodicShortestPathAnalysis(lattice(), config)
    search = interrupted._search

    def crash_in_period_2(graph, sources, period):
        if period == 2:
            raise KeyboardInterrupt
        return search(graph, sources, period)

    interrupted._search = crash_in_period_2
    with pytest.raises(KeyboardInterrupt):
        interrupted.run()
    assert checkpoint.exists()

    edited = PeriodicShortestPathAnalysis(lattice(shift=1e-6), config)
    assert edited._load_checkpoint(edited._checkpoint_params(float(N))) is None

    resumed = PeriodicShortestPathAnalysis(lattice(), config).run()
    fresh = PeriodicShortestPathAnalysis(lattice(), ShortestPathConfig(max_period=3)).run()
    assert np.array_equal(resumed.lengths, fresh.lengths)
    assert all(np.array_equal(a, b) for a, b in zip(resumed.paths, fresh.paths))
    assert not checkpoint.exists()  # deleted once the search completes


def test_saved_result(tmp_path):
    source = tmp_path / "network.data"
    source.write_text("")
    network = lattice()
    network = LammpsNetwork(network.atom_ids, network.coordinates, network.bond_ids, network.bond_types,
                            network.edges, network.box_bounds, source)
    result = PeriodicShortestPathAnalysis(network, ShortestPathConfig(max_period=3)).run()
    output = result.save(tmp_path / "network.x.sp.pkl")

    umask = os.umask(0)
    os.umask(umask)
    assert output.stat().st_mode & 0o777 == 0o666 & ~umask
    with open(output, "rb") as stream:
        loaded = pickle.load(stream)
    assert isinstance(loaded.source, str) and loaded.summary() == result.summary()
    with pytest.raises(FileExistsError):
        result.save(output)
    with pytest.raises(ValueError, match="must end in .pkl"):
        result.save(tmp_path / "network.x.sp.txt")
