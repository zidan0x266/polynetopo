"""Minimum cycle basis loop detection.

Uses igraph (much faster than networkx for this) to find the network's
minimum cycle basis, then counts how many crosslinkers sit on each loop -
a loop with exactly one crosslinker on it is a primary loop.

Reproducibility note: a minimum cycle basis isn't always unique - when
multiple equally-short cycles compete to cover the same part of the cycle
space, igraph's tie-breaking depends on the *order edges are added*, not on
how many unused/isolated vertex ids you padded the graph out to (tested
directly: on the ATRP example, padding the vertex count from 103149 - the
largest atom id actually referenced by a bond - up through 106250, one at a
time, never changed a single loop). So `n_vertices` is safe to infer; just
make sure the edges themselves are fed in a consistent order if you need to
reproduce a prior run's exact loop list byte-for-byte (`from_lammps` and
`from_polymer_network` both preserve the topology file's original line
order for this reason - reordering them, e.g. via `networkx`'s own edge
iteration, measurably changes a small number (13 of 2721 on the ATRP
example) of which loops get picked, though never the total count or the
primary-loop count). To make this a non-issue rather than something to
remember, both `LoopAnalyzer` constructors below funnel through the same
`PolymerNetwork` parse and its stored `.edges` - there is exactly one place
in this package that turns a topology file into an edge list, and both the
`networkx` and `igraph` views are built from that same list object, so they
can't disagree about node identity or edge order.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import igraph as ig

from .network import PolymerNetwork


@dataclass(slots=True)
class LoopBasis:
    """A network's minimum cycle basis, with a crosslinker count per loop."""

    loops: list[list[int]]
    crosslinker_counts: list[int]

    def as_legacy_tuple(self) -> list:
        """`[loops, counts]`, the shape earlier versions of this pipeline pickled."""
        return [self.loops, self.crosslinker_counts]

    def __len__(self) -> int:
        return len(self.loops)


class LoopAnalyzer:
    """Finds the minimum cycle basis of a network and its crosslinker loops."""

    def __init__(self, graph: ig.Graph) -> None:
        self.graph = graph

    @classmethod
    def from_lammps(cls, path: str | Path, n_vertices: Optional[int] = None) -> "LoopAnalyzer":
        """Build directly from a LAMMPS bond table.

        A convenience wrapper around `from_polymer_network`, for callers
        (like the standalone `loops` CLI subcommand) that don't already have
        a `PolymerNetwork` around. Parses the file through
        `PolymerNetwork.from_lammps` - the same parser `networkx` uses
        everywhere else in this package - rather than re-implementing its
        own, so there is exactly one LAMMPS-bond-table reader in the whole
        codebase and no way for its edge order to drift out of sync between
        the `networkx` and `igraph` views of the same topology.
        """
        return cls.from_polymer_network(PolymerNetwork.from_lammps(path), n_vertices)

    @classmethod
    def from_polymer_network(
        cls, network: PolymerNetwork, n_vertices: Optional[int] = None
    ) -> "LoopAnalyzer":
        """Build from an already-loaded `PolymerNetwork`, without re-reading its file.

        Uses `network.edges` (the original file order) rather than
        `network.graph.edges()` (`networkx`'s own iteration order, which
        generally differs once a graph has any branching) - see the module
        docstring for why edge order, not vertex count, is what matters here.
        """
        if n_vertices is None:
            n_vertices = max(network.graph.nodes())

        graph = ig.Graph(directed=False)
        graph.add_vertices(n_vertices + 1)
        graph.add_edges(network.edges)
        return cls(graph)

    def minimum_cycle_basis(self, crosslinker_degree: int = 3) -> LoopBasis:
        """The minimum cycle basis, with per-loop crosslinker counts.

        A loop's crosslinker count is the number of its nodes with degree
        >= `crosslinker_degree`; a count of 1 marks a primary loop.
        """
        cycles = self.graph.minimum_cycle_basis()

        loops: list[list[int]] = []
        for cycle in cycles:
            loop: list[int] = []
            for edge_idx in cycle:
                edge = self.graph.es[edge_idx]
                for node in (edge.source, edge.target):
                    if node not in loop:
                        loop.append(node)
            loops.append(loop)

        counts = [
            sum(1 for node in loop if self.graph.vs[node].degree() >= crosslinker_degree)
            for loop in loops
        ]
        return LoopBasis(loops=loops, crosslinker_counts=counts)
