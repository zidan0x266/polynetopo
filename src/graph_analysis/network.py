"""Loading polymer network topologies into graphs."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterator, Optional

import networkx as nx
import numpy as np


class PolymerNetwork:
    """An undirected polymer network topology.

    Thin wrapper around a :class:`networkx.Graph` that knows how to build
    itself from the edge-list formats used in this pipeline, and exposes the
    small set of graph queries the rest of the pipeline needs.

    `from_lammps` (a full LAMMPS data file's Bonds section) is the primary,
    expected format. `from_rmc` and `from_bond_list` are legacy formats kept
    for older topology files that were pre-trimmed to a bare bond list
    before this package could read a full data file directly.

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
    def from_lammps(cls, path: str | Path) -> "PolymerNetwork":
        """Load a graph from a full LAMMPS data file (the primary format).

        Reads only the Bonds section of a `write_data` output (with its
        Masses/Atoms/Velocities/Bonds/... sections) - the rest (coordinates,
        masses, etc.) isn't needed for topology analysis. If you have an
        already-trimmed bare bond table instead, use `from_bond_list`.
        """
        return cls._from_edges(_read_lammps_data_bonds_section(path), path)

    @classmethod
    def from_bond_list(cls, path: str | Path) -> "PolymerNetwork":
        """Legacy: load a graph from a bare bond table.

        Each line is `bond-id bond-type atom_i atom_j` and nothing else -
        what manually trimming a LAMMPS data file's Bonds section out to its
        own file produces. Kept for topology files saved this way before
        `from_lammps` could read a full data file directly; also a
        reasonable starting point to adapt for a bond list produced by a
        different tool (e.g. a reactive Monte Carlo / 3DBFM code) that isn't
        LAMMPS-shaped, if its column layout differs from this one.
        """
        return cls._from_edges(_read_bond_table(path), path)

    @classmethod
    def from_rmc(cls, path: str | Path) -> "PolymerNetwork":
        """Legacy: load a graph from an RMC-style whitespace edge list.

        Each line holds an even number of atom ids; consecutive pairs are
        edges. A bond may be listed in both directions or in one: of a pair
        listed both ways only the direction with the lower ID first is kept.
        """
        pairs: list[tuple[int, int]] = []
        with open(path) as data_file:
            for line in data_file:
                ids = [int(value) for value in line.split()]
                pairs.extend(zip(ids[::2], ids[1::2]))
        listed = set(pairs)
        edges = [(a, b) for a, b in pairs if a <= b or (b, a) not in listed]
        return cls._from_edges(edges, path)

    @classmethod
    def _from_edges(cls, edges: list[tuple[int, int]], path: str | Path) -> "PolymerNetwork":
        _require_simple(edges, path)
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


def _require_simple(edges: list[tuple[int, int]], path: str | Path) -> None:
    """Reject a bond from an atom to itself or an atom pair bonded more than once."""
    seen: set[tuple[int, int]] = set()
    repeated: list[tuple[int, int]] = []
    for a, b in edges:
        if a == b:
            raise ValueError(f"{path}: atom {a} is bonded to itself")
        pair = (a, b) if a < b else (b, a)
        if pair in seen:
            repeated.append(pair)
        seen.add(pair)
    if repeated:
        a, b = repeated[0]
        raise ValueError(
            f"{path}: {len(repeated)} atom pairs are bonded more than once (first: {a} and {b}); "
            "a strand lists atoms, so it cannot tell parallel bonds apart"
        )


def _read_bond_table(path: str | Path) -> list[tuple[int, int]]:
    data = np.loadtxt(path, ndmin=2)
    return [(int(atom_i), int(atom_j)) for atom_i, atom_j in data[:, (2, 3)]]


def _read_lammps_data_bonds_section(path: str | Path) -> list[tuple[int, int]]:
    """Extract just the Bonds section of a full LAMMPS data file.

    Sections are separated by blank lines, each headed by its name (e.g.
    "Masses", "Atoms # full", "Bonds") on its own line. Each bond data line
    is `bond-id bond-type atom_i atom_j`.
    """
    edges: list[tuple[int, int]] = []
    in_bonds = False
    with open(path) as f:
        for line in f:
            stripped = line.split("#", 1)[0].strip()
            if not in_bonds:
                if stripped == "Bonds":
                    in_bonds = True
                continue
            if not stripped:
                if edges:
                    break  # blank line after data rows: end of the section
                continue  # the blank line right after the "Bonds" header
            _, _, atom_i, atom_j = stripped.split()
            edges.append((int(atom_i), int(atom_j)))

    if not edges:
        raise ValueError(
            f"No 'Bonds' section found in {path}. If this is already a bare bond "
            "table (bond-id bond-type atom_i atom_j, nothing else), use "
            "PolymerNetwork.from_bond_list() instead."
        )
    return edges
