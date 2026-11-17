"""GEBC/OBC/DBC of the in-house igraph build; skipped with stock igraph."""

from __future__ import annotations

import igraph as ig
import numpy as np
import pytest

from graph_analysis.lammps import LammpsNetwork

from test_shortest_path import N, lattice

pytestmark = pytest.mark.skipif(not hasattr(ig.Graph, "edge_betweenness_spatial"),
                                reason="needs the in-house igraph build in vendor/")


def lattice_with_chain_end() -> LammpsNetwork:
    """The lattice of test_shortest_path plus one bead bonded to bead 0 only."""
    base = lattice()
    end = base.n_atoms
    coordinates = np.vstack([base.coordinates, base.coordinates[0] + (0.3, 0.2, 0.1)])
    edges = np.vstack([base.edges, (0, end)])
    return LammpsNetwork(np.arange(1, end + 2), coordinates, np.arange(1, len(edges) + 1), np.ones(len(edges), int),
                         edges, base.box_bounds)


def test_chain_end_bond_has_closed_forms():
    from graph_analysis.betweenness import AnalysisConfig, SpatialBetweennessAnalysis

    network = lattice_with_chain_end()
    result = SpatialBetweennessAnalysis(network, AnalysisConfig(direction="x", validate_stock=True)).run()
    assert result.validation.allclose  # native GEBC equals stock igraph's

    # every shortest path from the chain end runs through its bond: one pair per other bead,
    # weighted by the minimum-image separation (DBC) and its alignment with x (OBC)
    end, bond = network.n_atoms - 1, network.n_bonds - 1
    d = np.delete(network.coordinates - network.coordinates[end], end, axis=0)
    d -= N * np.round(d / N)
    r = np.linalg.norm(d, axis=1)
    assert result.raw_gebc[bond] == network.n_atoms - 1
    assert np.isclose(result.raw_dbc[bond], r.sum(), rtol=1e-12)
    assert np.isclose(result.raw_obc[bond], (np.abs(d[:, 0]) / r).sum(), rtol=1e-12)
