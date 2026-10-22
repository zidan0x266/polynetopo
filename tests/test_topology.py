"""Strand/dangling search and the topology stages of the command line."""

from __future__ import annotations

import pickle

import networkx as nx
import pytest

from graph_analysis import cli
from graph_analysis.crosslink import CrosslinkAction, CrosslinkPathAnalyzer
from graph_analysis.network import PolymerNetwork


def network(edges):
    graph = nx.Graph()
    graph.add_edges_from(edges)
    return PolymerNetwork(graph, edges)


# Crosslinkers 1 and 2 are joined by three separate chains - a direct bond, one atom (3) in the
# middle, and two atoms (4, 5) - and crosslinker 1 also carries a dangling end (6). An earlier
# pairwise search missed the third chain of any pair joined by three or more.
PARALLEL = [(1, 2), (1, 3), (3, 2), (1, 4), (4, 5), (5, 2), (1, 6)]


def test_every_parallel_strand_is_found():
    strands = CrosslinkPathAnalyzer(CrosslinkAction.CROSSLINKERS).run(network(PARALLEL))
    assert sorted(strands) == [[1, 2], [1, 3, 2], [1, 4, 5, 2]]


def test_dangling_chain_runs_from_a_crosslinker_to_a_chain_end():
    assert CrosslinkPathAnalyzer(CrosslinkAction.DANGLINGS).run(network(PARALLEL)) == [[1, 6]]


def test_paths_run_from_the_lower_to_the_higher_atom_id():
    strands = CrosslinkPathAnalyzer(CrosslinkAction.CROSSLINKERS).run(network(PARALLEL + [(2, 0), (0, 7)]))
    assert all(path[0] < path[-1] for path in strands)


def load(path):
    with open(path, "rb") as stream:
        return pickle.load(stream)


def test_full_stage_writes_the_files_that_scanlan_case_reads(random_file, tmp_path):
    results = tmp_path / "results.pkl"
    cli.main(["--function", "topology", "--input-file", str(random_file), "--output-file", str(results), "--raw"])
    for suffix in ("strands", "danglings", "loops", "raw"):
        assert (tmp_path / f"results_{suffix}.pkl").exists()

    # each single stage gives the same file as the full run wrote
    cli.main(["--function", "topology", "--stage", "strands", "--input-file", str(random_file),
              "--output-file", str(tmp_path / "s.pkl")])
    cli.main(["--function", "topology", "--stage", "danglings", "--input-file", str(random_file),
              "--output-file", str(tmp_path / "d.pkl")])
    cli.main(["--function", "topology", "--stage", "loops", "--input-file", str(random_file),
              "--output-file", str(tmp_path / "l.pkl")])
    assert load(tmp_path / "s.pkl") == load(tmp_path / "results_strands.pkl")
    assert load(tmp_path / "d.pkl") == load(tmp_path / "results_danglings.pkl")
    assert load(tmp_path / "l.pkl") == load(tmp_path / "results_loops.pkl")

    # the refinement alone, from those files, reproduces the full run's refined result
    cli.main(["--function", "topology", "--stage", "scanlan-case", "--input-file", str(random_file),
              "--strands-file", str(tmp_path / "results_strands.pkl"),
              "--danglings-file", str(tmp_path / "results_danglings.pkl"),
              "--loops-file", str(tmp_path / "results_loops.pkl"),
              "--output-file", str(tmp_path / "again.pkl")])
    assert load(tmp_path / "again.pkl") == load(results)


def test_a_single_stage_never_runs_the_refinement(random_file, tmp_path):
    cli.main(["--function", "topology", "--stage", "strands", "--input-file", str(random_file),
              "--output-file", str(tmp_path / "s.pkl")])
    assert [p.name for p in tmp_path.glob("*.pkl")] == ["s.pkl"]


def test_scanlan_case_stage_needs_its_three_files(random_file, tmp_path):
    with pytest.raises(SystemExit):
        cli.main(["--function", "topology", "--stage", "scanlan-case", "--input-file", str(random_file),
                  "--output-file", str(tmp_path / "x.pkl")])


def test_removed_options_are_rejected(random_file):
    for argv in (["--function", "chains"], ["--function", "shortest-path", "--image-max", "8"]):
        with pytest.raises(SystemExit):
            cli.main(argv + ["--input-file", str(random_file)])
