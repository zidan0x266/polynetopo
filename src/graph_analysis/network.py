"""Loading polymer network topologies into graphs."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterator, Optional

import networkx as nx
import numpy as np


class PolymerNetwork:
    """An undirected polymer network topology.

    Thin wrapper around a :class:`networkx.Graph` that knows how to build
    itself from the edge-list formats used in this pipeline (RMC-style text
    files and LAMMPS bond tables), and exposes the small set of graph
    queries the rest of the pipeline needs.

    Also keeps the original file-order edge list separately from the graph:
    `networkx`'s own edge iteration order doesn't generally preserve
    insertion order once a graph has any branching, and igraph's minimum
    cycle basis (see `loops.py`) can pick a different - equally valid, but
    different - basis depending on the order edges are fed to it. Loaders
    that need a specific prior result reproduced should use `.edges`, not
    `.graph.edges()`.
    """

    __slots__ = ("graph", "edges")

    def __init__(self, graph: nx.Graph, edges: Optional[list[tuple[int, int]]] = None) -> None:
        self.graph = graph
        self.edges = edges if edges is not None else list(graph.edges())

    @classmethod
    def from_rmc(cls, path: str | Path) -> "PolymerNetwork":
        """Load a graph from an RMC-style whitespace edge list.

        Each line holds an even number of atom ids; consecutive pairs are
        edges. Only one direction of each unordered pair is kept.
        """
        edges: list[tuple[int, int]] = []
        with open(path) as data_file:
            for line in data_file:
                ids = line.split()
                atom_a, atom_b = ids[::2], ids[1::2]
                edges.extend(
                    (int(a), int(b)) for a, b in zip(atom_a, atom_b) if int(a) < int(b)
                )
        graph = nx.Graph()
        graph.add_edges_from(edges)
        return cls(graph, edges)

    @classmethod
    def from_lammps(cls, path: str | Path) -> "PolymerNetwork":
        """Load a graph from a LAMMPS bond table (atom ids in columns 3 and 4)."""
        data = np.loadtxt(path)
        edges = [(int(atom_i), int(atom_j)) for atom_i, atom_j in data[:, (2, 3)]]
        graph = nx.Graph()
        graph.add_edges_from(edges)
        return cls(graph, edges)

    def nodes_with_degree(self, predicate: Callable[[int], bool]) -> list[int]:
        """Nodes whose degree satisfies `predicate`."""
        return [node for node, degree in self.graph.degree() if predicate(degree)]

    def connected_components(self) -> Iterator["PolymerNetwork"]:
        """Yield one `PolymerNetwork` per connected component."""
        for component in nx.connected_components(self.graph):
            yield PolymerNetwork(self.graph.subgraph(component).copy())

    def __len__(self) -> int:
        return self.graph.number_of_nodes()

    def __repr__(self) -> str:
        return f"PolymerNetwork(nodes={self.graph.number_of_nodes()}, edges={self.graph.number_of_edges()})"
