"""The shipped examples: the refinement of example 1 from its stage files, and the results of examples 2 and 3."""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pytest

from graph_analysis import LoopBasis, NetworkStats, PolymerNetwork, ScanlanCaseAnalyzer

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def load(path: Path):
    with open(path, "rb") as stream:
        return pickle.load(stream)


def test_example_1_refinement():
    folder = EXAMPLES / "example_1"
    network = PolymerNetwork.from_lammps(EXAMPLES / "data" / "atrp_s1_sys1_0.data")
    raw = NetworkStats.from_loop_basis(load(folder / "results_strands.pkl"), load(folder / "results_danglings.pkl"),
                                       LoopBasis(*load(folder / "results_loops.pkl")))
    refined = ScanlanCaseAnalyzer(network.graph, raw).run().refined
    counts = (refined.effective_crosslinker_d4, refined.effective_crosslinker_d3, len(refined.strands), len(refined.danglings))
    assert counts == (2006, 396, 4606, 487)  # as in TUTORIAL.md
    shipped = load(folder / "results.pkl")
    assert shipped.strands == refined.strands and shipped.danglings == refined.danglings
    assert shipped.loop_orders == refined.loop_orders


def test_example_2_betweenness():
    result = load(EXAMPLES / "example_2" / "atrp_s1_sys1_0.x.spatial.pkl")
    summary = json.loads((EXAMPLES / "example_2" / "summary.json").read_text())
    assert result.summary() == {key: value for key, value in summary.items() if key != "output"}
    assert result.raw_gebc.shape == (105869,)
    assert int(result.network.bond_ids[np.argmax(result.raw_gebc)]) == 69986  # as in TUTORIAL.md


def test_example_3_shortest_paths():
    result = load(EXAMPLES / "example_3" / "atrp_s1_sys1_0.x.sp.pkl")
    summary = json.loads((EXAMPLES / "example_3" / "summary.json").read_text())
    assert result.summary() == {key: value for key, value in summary.items() if key != "output"}
    assert isinstance(result.source, str)


@pytest.mark.parametrize("example", ["example_2/atrp_s1_sys1_0.x.spatial.pkl", "example_3/atrp_s1_sys1_0.x.sp.pkl"])
def test_pickles_do_not_reference_pathlib(example):
    assert b"pathlib" not in (EXAMPLES / example).read_bytes()  # loadable under any Python version
