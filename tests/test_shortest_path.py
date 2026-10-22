"""Periodic shortest paths (graph_analysis.shortest_path)."""

from __future__ import annotations

import pickle

import numpy as np
import pytest

from graph_analysis.shortest_path import PeriodicShortestPathAnalysis, ShortestPathConfig


def run(path, **settings):
    return PeriodicShortestPathAnalysis.from_lammps(path, ShortestPathConfig(**settings)).run()


def test_ring_has_a_single_shortest_path(ring_file):
    result = run(ring_file, max_period=4)
    # Every atom finds the same loop of 12 bonds; Case 1 keeps one copy. At p = 2 the shortest
    # route goes round twice and revisits atoms, so it is invalid and the search stops.
    assert result.lengths.tolist() == [12.0]
    assert result.min_length == 12.0
    assert result.min_length_by_period.tolist() == [[1.0, 12.0]]
    assert len(result.paths) == 1
    path = result.paths[0]
    assert path[0, 0] == path[-1, 0]  # starts and ends on the same atom
    assert path[0, 1] == 0 and path[-1, 1] == 1  # image 0 to image p = 1
    assert len(path) - 1 == 12


def test_keeping_redundant_paths_keeps_every_start_atom(ring_file):
    result = run(ring_file, max_period=1, remove_redundant=False)
    assert len(result.paths) == 12
    assert sorted(int(p[0, 0]) for p in result.paths) == list(range(1, 13))


def test_nothing_wraps_raises(chain_file):
    with pytest.raises(ValueError, match="No shortest path"):
        run(chain_file)


def test_repeated_bond_raises(ring_file):
    text = ring_file.read_text().rstrip("\n") + "\n13 1 2 1\n"  # atoms 1 and 2 are already bonded
    ring_file.write_text(text.replace("12 bonds", "13 bonds"))
    with pytest.raises(ValueError, match="bonded more than once"):
        PeriodicShortestPathAnalysis.from_lammps(ring_file)


def test_config_validation():
    for bad in ({"max_period": 0}, {"buffer_images": -1}, {"slice_max": 0}, {"ncpus": 0}):
        with pytest.raises(ValueError):
            ShortestPathConfig(**bad)
    assert ShortestPathConfig().max_period == 8  # the paper's M
    assert ShortestPathConfig().min_sp_only is False  # Case 3 is optional
    assert ShortestPathConfig(direction="z").direction == 2


# Counts and minimum D^g checked against the reference code (the Jackson lab's shortest_path_thermoset,
# as run on the cluster) on the `random_file` network, for M = 4.
@pytest.mark.parametrize(
    "direction, buffer, paths, case_3_paths, minimum",
    [
        ("x", 0, 169, 117, 12.0),
        ("x", 1, 171, 118, 12.0),
        ("y", 0, 209, 131, 8.0),
        ("y", 1, 209, 132, 8.0),
        ("z", 0, 134, 114, 10.0),
        ("z", 1, 131, 114, 10.0),
    ],
)
def test_matches_reference_on_random_network(random_file, direction, buffer, paths, case_3_paths, minimum):
    settings = dict(direction=direction, buffer_images=buffer, max_period=4)
    result = run(random_file, **settings)
    assert len(result.paths) == paths
    assert result.min_length == minimum
    case_3 = run(random_file, **settings, min_sp_only=True)
    assert len(case_3.paths) == case_3_paths
    assert case_3.min_length == minimum  # Case 3 never changes the first (shortest) SP


def test_min_sp_only_keeps_the_shortest_path_of_each_atom(random_file):
    everything = run(random_file, max_period=4, remove_redundant=False)
    shortest = run(random_file, max_period=4, remove_redundant=False, min_sp_only=True)
    starts = [int(p[0, 0]) for p in shortest.paths]
    assert starts == sorted(set(starts))  # one path per atom, in atom order
    best = {}
    for length, path in zip(everything.lengths, everything.paths):
        atom = int(path[0, 0])
        best[atom] = min(best.get(atom, np.inf), length)
    assert shortest.lengths.tolist() == [best[atom] for atom in starts]


def test_max_period_limits_the_search(random_file):
    one = run(random_file, max_period=1)
    assert one.min_length_by_period[:, 0].tolist() == [1.0]
    four = run(random_file, max_period=4)
    assert max(four.min_length_by_period[:, 0]) <= 4
    # the periods searched first do not depend on how far the search is allowed to go
    assert four.min_length_by_period[0].tolist() == one.min_length_by_period[0].tolist()


def test_paths_start_in_image_zero_and_end_on_the_same_atom(random_file):
    for path in run(random_file, max_period=3).paths:
        assert path[0, 0] == path[-1, 0]
        assert path[0, 1] == 0
        assert path[-1, 1] >= 1


def test_workers_give_the_same_result_as_one_process(random_file):
    serial = run(random_file, max_period=3)
    parallel = run(random_file, max_period=3, ncpus=2)
    assert np.array_equal(serial.lengths, parallel.lengths)
    assert all(np.array_equal(a, b) for a, b in zip(serial.paths, parallel.paths))


def test_slice_only_uses_nearby_start_atoms(random_file):
    full = run(random_file, max_period=2)
    sliced = run(random_file, max_period=2, slice_max=3.0)
    assert len(sliced.paths) < len(full.paths)
    assert sliced.min_length >= full.min_length


def test_save_and_load_round_trip(ring_file, tmp_path):
    result = run(ring_file, max_period=2)
    output = result.save(tmp_path / "ring.x.sp.pkl")
    with open(output, "rb") as stream:
        loaded = pickle.load(stream)
    assert np.array_equal(loaded.lengths, result.lengths)
    assert loaded.summary()["max_period"] == 2
    with pytest.raises(FileExistsError):
        result.save(output)
    result.save(output, overwrite=True)
    with pytest.raises(ValueError):
        result.save(tmp_path / "wrong_suffix.txt")
