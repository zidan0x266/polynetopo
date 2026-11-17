"""Scanlan-Case refinement of strand/dangling/loop statistics.

Per the Scanlan-Case criterion (Zhang et al. 2025, *Macromolecules*,
"Topology and Mechanical Properties of Polymer Networks Formed under Free
Radical and Atom Transfer Radical Polymerizations"), a crosslinker is
elastically effective only if at least three of its arms lead to the network
independently of one another. For a single defect this means: a crosslinker
of functionality f=4 that carries a primary loop or two dangling chains is a
linear unit, and one that carries one dangling chain is an f=3 junction; a
crosslinker of f=3 that carries a primary loop or a dangling chain is a
linear unit, folded into whichever strand runs through it.

The criterion is applied to the whole network at once, so that defects that
build on each other are included as well (a strand that leads only to a dead
end is itself dangling, and the crosslinker at its other end loses an arm):
the load-bearing part of the network is its largest biconnected component,
i.e., the network without everything that hangs on it through a single strand
or a single crosslinker - dangling chains, primary loops, and clusters of
crosslinkers that lead only to a dead end. The effective functionality of a
crosslinker is its number of bonds inside that component. This module
reclassifies the raw strands/danglings/loops accordingly and recomputes the
effective crosslinker density.
"""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from dataclasses import dataclass, field

import networkx as nx

from .loops import LoopBasis

logger = logging.getLogger(__name__)

Strand = list[int]


def merge_strands(node: int, strand_a: Strand, strand_b: Strand) -> Strand:
    """Splice two strands that share an endpoint at `node` into one strand."""
    a = strand_a if strand_a[-1] == node else strand_a[::-1]
    b = strand_b if strand_b[0] == node else strand_b[::-1]
    return a + b[1:]


def load_bearing_core(graph: nx.Graph) -> set[int]:
    """Atoms of the largest biconnected component, or none if it has no crosslinker of degree >= 3 in it."""
    core = max(nx.biconnected_components(graph), key=len, default=set())
    if not any(sum(1 for n in graph.adj[node] if n in core) >= 3 for node in core):
        return set()
    return core


@dataclass(slots=True)
class NetworkStats:
    """Strand/dangling/loop statistics for one snapshot of a network (raw or refined)."""

    strands: list[Strand]
    danglings: list[Strand]
    loops: list[list[int]]
    loop_orders: list[int]

    @classmethod
    def from_loop_basis(
        cls, strands: list[Strand], danglings: list[Strand], loops: LoopBasis
    ) -> "NetworkStats":
        return cls(
            strands=strands,
            danglings=danglings,
            loops=loops.loops,
            loop_orders=loops.crosslinker_counts,
        )


@dataclass(slots=True)
class RefinedStats(NetworkStats):
    """`NetworkStats` after Scanlan-Case refinement, plus the crosslinker density it produced.

    `effective_crosslinker_d4/d3` count the crosslinkers of effective functionality
    4 and 3, `all_crosslinker_d4/d3` the atoms of degree 4 and 3, and
    `correct_crosslinker_d4/d3` = all - effective, the correction to each count;
    `correct_crosslinker_d3` is negative when more f=4 crosslinkers become f=3
    junctions (`switched_crosslinker_d4`) than f=3 crosslinkers become linear units.
    `effective_functionality` maps every effective crosslinker to its effective
    functionality, including any above 4.
    """

    effective_crosslinker_d4: int = 0
    effective_crosslinker_d3: int = 0
    all_crosslinker_d4: int = 0
    all_crosslinker_d3: int = 0
    correct_crosslinker_d4: int = 0
    correct_crosslinker_d3: int = 0
    switched_crosslinker_d4: int = 0
    effective_functionality: dict[int, int] = field(default_factory=dict)


@dataclass(slots=True)
class ScanlanCaseResult:
    """Network statistics before and after applying the Scanlan-Case criterion."""

    refined: RefinedStats
    raw: NetworkStats

    def to_legacy(self) -> list:
        """The nested-list form `[[refined...], [raw...], [counts...]]`."""
        return [
            [self.refined.strands, self.refined.danglings, self.refined.loops, self.refined.loop_orders],
            [self.raw.strands, self.raw.danglings, self.raw.loops, self.raw.loop_orders],
            [
                self.refined.effective_crosslinker_d4,
                self.refined.effective_crosslinker_d3,
                self.refined.all_crosslinker_d4,
                self.refined.all_crosslinker_d3,
                self.refined.correct_crosslinker_d4,
                self.refined.correct_crosslinker_d3,
                self.refined.switched_crosslinker_d4,
            ],
        ]


class ScanlanCaseAnalyzer:
    """Applies the Scanlan-Case criterion to a network's strands/danglings/loops.

    The refined strands are the raw strands inside the load-bearing component,
    with the strands that meet at a crosslinker of effective functionality 2 (a
    linear unit) spliced into one; every other raw strand is added to the
    dangling chains, so each raw strand ends up in exactly one of the two
    lists. A loop keeps order 1 if it is a primary loop and otherwise counts
    the effective crosslinkers on it.
    """

    def __init__(self, graph: nx.Graph, raw: NetworkStats) -> None:
        self.graph = graph
        self.raw = raw

    def run(self) -> ScanlanCaseResult:
        graph, raw = self.graph, self.raw

        # Load-bearing component and the effective functionality inside it.
        core = load_bearing_core(graph)
        core_degree = {node: sum(1 for n in graph.adj[node] if n in core) for node in core}
        functionality = {node: degree for node, degree in core_degree.items() if degree >= 3}
        linear = {node for node, degree in core_degree.items() if degree == 2 and graph.degree(node) >= 3}

        # Strands inside the component are kept, spliced at its linear units;
        # the others lead to a dead end and become dangling.
        kept, removed = [], []
        for strand in raw.strands:
            (kept if all(node in core for node in strand) else removed).append(strand)
        refined_strands = _splice(kept, linear)
        refined_danglings = raw.danglings + removed
        logger.info(
            "Load-bearing component: %d atoms, %d effective crosslinkers; %d strands, %d moved to the danglings",
            len(core), len(functionality), len(refined_strands), len(removed),
        )

        refined_loop_orders = [
            order if order == 1 else sum(1 for node in loop if node in functionality)
            for order, loop in zip(raw.loop_orders, raw.loops)
        ]

        degree_counts = Counter(degree for _, degree in graph.degree)
        effective_counts = Counter(functionality.values())
        all_d4, all_d3 = degree_counts[4], degree_counts[3]
        effective_d4, effective_d3 = effective_counts[4], effective_counts[3]

        return ScanlanCaseResult(
            refined=RefinedStats(
                strands=refined_strands,
                danglings=refined_danglings,
                loops=raw.loops,
                loop_orders=refined_loop_orders,
                effective_crosslinker_d4=effective_d4,
                effective_crosslinker_d3=effective_d3,
                all_crosslinker_d4=all_d4,
                all_crosslinker_d3=all_d3,
                correct_crosslinker_d4=all_d4 - effective_d4,
                correct_crosslinker_d3=all_d3 - effective_d3,
                switched_crosslinker_d4=sum(
                    1 for node, f in functionality.items() if f == 3 and graph.degree(node) == 4
                ),
                effective_functionality=functionality,
            ),
            raw=raw,
        )


def _splice(strands: list[Strand], linear: set[int]) -> list[Strand]:
    """Join the strands that meet at each node of `linear` (exactly two do), lower end ID first."""
    at_node: dict[int, list[int]] = defaultdict(list)
    for index, strand in enumerate(strands):
        for end in (strand[0], strand[-1]):
            if end in linear:
                at_node[end].append(index)

    used = [False] * len(strands)
    spliced: list[Strand] = []
    for index, strand in enumerate(strands):
        if used[index]:
            continue
        used[index] = True
        path = list(strand)
        for _ in range(2):  # extend the last end, then turn round and extend the other one
            while path[-1] in linear:
                node = path[-1]
                following = next((i for i in at_node[node] if not used[i]), None)
                if following is None:
                    break
                used[following] = True
                path = merge_strands(node, path, strands[following])
            path.reverse()
        spliced.append(path if path[0] <= path[-1] else path[::-1])
    return spliced
