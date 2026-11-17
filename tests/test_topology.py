"""Topology: the Scanlan-Case refinement on hand-built networks, and the checks of the topology input."""

from __future__ import annotations

import pytest

from graph_analysis import LoopAnalyzer, NetworkAnalysisPipeline, PolymerNetwork
from graph_analysis.cli import main

# The load-bearing part of every case below: four crosslinkers joined pairwise by six two-bond strands.
CORE = [[10, 11, 20], [10, 12, 30], [10, 13, 40], [20, 14, 30], [20, 15, 40], [30, 16, 40]]
K4 = [bond for strand in CORE for bond in zip(strand, strand[1:])]

# Defects added to the core. Each hangs on it through a single strand or a single crosslinker, so the refined
# network is the core itself, with node 10 left at f=3.
DEFECTS = {
    "primary loop and dangling chain on one crosslinker": [(1, 2), (2, 3), (3, 1), (1, 4), (1, 5), (5, 10)],
    "f=4 crosslinker with three dangling chains": [(1, 2), (1, 3), (1, 4), (1, 5), (5, 10)],
    "two crosslinkers that lead only to dangling chains": [(1, 2), (1, 10), (1, 3), (3, 4), (4, 5), (4, 6)],
    "cluster hanging on one strand": [(1, 2), (1, 3), (3, 2), (1, 4), (4, 2), (1, 5), (5, 10)],
    "cluster hanging on one crosslinker": [(1, 2), (1, 3), (3, 2), (1, 4), (4, 2), (1, 10), (2, 10)],
}


def bond_list(tmp_path, name, edges):
    path = tmp_path / f"{name}.bds"
    path.write_text("".join(f"{i + 1} 1 {a} {b}\n" for i, (a, b) in enumerate(edges)))
    return path


def refine(tmp_path, edges):
    network = PolymerNetwork.from_bond_list(bond_list(tmp_path, "network", edges))
    pipeline = NetworkAnalysisPipeline(network)
    return pipeline.raw_stats(), pipeline.run().refined


def bonds(paths):
    return sum(len(path) - 1 for path in paths)


@pytest.mark.parametrize("defect", list(DEFECTS))
def test_defects_leave_the_load_bearing_core(tmp_path, defect):
    raw, refined = refine(tmp_path, K4 + DEFECTS[defect])
    assert sorted(refined.strands) == CORE
    assert (refined.effective_crosslinker_d4, refined.effective_crosslinker_d3) == (0, 4)
    assert refined.effective_functionality == {10: 3, 20: 3, 30: 3, 40: 3}
    # every raw strand ends up in exactly one of the two lists
    assert bonds(refined.strands) + bonds(refined.danglings) == bonds(raw.strands) + bonds(raw.danglings)


def test_crosslinker_of_functionality_five(tmp_path):
    # node 1 carries a primary loop and three strands to K4: f=5 -> f=3, and 10, 20, 30 become f=4
    raw, refined = refine(tmp_path, K4 + [(1, 2), (2, 3), (3, 1), (1, 4), (4, 10), (1, 5), (5, 20), (1, 6), (6, 30)])
    assert len(refined.strands) == 9
    assert refined.effective_functionality == {1: 3, 10: 4, 20: 4, 30: 4, 40: 3}
    assert (refined.all_crosslinker_d4, refined.correct_crosslinker_d4, refined.correct_crosslinker_d3) == (3, 0, -1)


def test_strands_are_spliced_at_linear_units(tmp_path):
    # node 1 joins 10 and 20 and carries a dangling chain: a linear unit, so 10-1-20 is one strand
    raw, refined = refine(tmp_path, K4 + [(10, 1), (1, 20), (1, 2)])
    assert [10, 1, 20] in refined.strands
    assert 1 not in refined.effective_functionality
    assert refined.effective_functionality[10] == 4 and refined.effective_functionality[20] == 4


@pytest.mark.parametrize("extra, message", [([(10, 11)], "bonded more than once"), ([(10, 10)], "bonded to itself")])
def test_repeated_and_self_bonds_are_rejected(tmp_path, extra, message):
    with pytest.raises(ValueError, match=message):
        PolymerNetwork.from_bond_list(bond_list(tmp_path, "bad", K4 + extra))


def test_one_line_bond_table(tmp_path):
    path = tmp_path / "one.bds"
    path.write_text("1 1 3 7\n")
    assert PolymerNetwork.from_bond_list(path).edges == [(3, 7)]


def test_rmc_bonds_listed_once_or_both_ways(tmp_path):
    once, both = tmp_path / "once.rmc", tmp_path / "both.rmc"
    once.write_text("20 10 30 10\n40 10\n")
    both.write_text("10 20 10 30\n20 10 30 10\n")
    assert PolymerNetwork.from_rmc(once).edges == [(20, 10), (30, 10), (40, 10)]
    assert PolymerNetwork.from_rmc(both).edges == [(10, 20), (10, 30)]


def test_nvertices_is_the_largest_atom_id(tmp_path):
    network = PolymerNetwork.from_bond_list(bond_list(tmp_path, "k4", K4))
    with pytest.raises(ValueError, match="largest atom ID"):
        LoopAnalyzer.from_polymer_network(network, n_vertices=4)
    assert len(LoopAnalyzer.from_polymer_network(network).minimum_cycle_basis()) == 3


def test_topology_keeps_existing_output_without_overwrite(tmp_path):
    argv = ["--function", "topology", "--format", "bond-list", "--input-file", str(bond_list(tmp_path, "k4", K4)),
            "--output-file", str(tmp_path / "results.pkl")]
    main(argv)
    with pytest.raises(FileExistsError):
        main(argv)
    main(argv + ["--overwrite"])
