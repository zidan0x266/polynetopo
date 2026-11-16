"""Scanlan-Case refinement of strand/dangling/loop statistics.

Per the Scanlan-Case criterion (Zhang et al. 2025, *Macromolecules*,
"Topology and Mechanical Properties of Polymer Networks Formed under Free
Radical and Atom Transfer Radical Polymerizations"): a crosslinker of
functionality f=4 that carries a primary loop or two dangling chains isn't
actually load-bearing, and is reclassified as a linear unit; a crosslinker
of f=3 that carries a primary loop or a dangling chain is reclassified as a
linear unit, folded into whichever dangling chain or strand runs through it.
This module reclassifies the raw strands/danglings/loops accordingly and
recomputes the effective crosslinker density.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import networkx as nx

from .loops import LoopBasis

logger = logging.getLogger(__name__)

Strand = list[int]


def merge_strands(node: int, strand_a: Strand, strand_b: Strand) -> Strand:
    """Splice two strands that share an endpoint at `node` into one strand."""
    a = strand_a if strand_a[-1] == node else strand_a[::-1]
    b = strand_b if strand_b[0] == node else strand_b[::-1]
    return a + b[1:]


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
    """`NetworkStats` after Scanlan-Case refinement, plus the crosslinker density it produced."""

    effective_crosslinker_d4: int = 0
    effective_crosslinker_d3: int = 0
    all_crosslinker_d4: int = 0
    all_crosslinker_d3: int = 0
    correct_crosslinker_d4: int = 0
    correct_crosslinker_d3: int = 0
    switched_crosslinker_d4: int = 0


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
    """Applies the Scanlan-Case criterion to a network's strands/danglings/loops."""

    def __init__(self, graph: nx.Graph, raw: NetworkStats) -> None:
        self.graph = graph
        self.raw = raw

    def run(self) -> ScanlanCaseResult:
        graph = self.graph
        strands, danglings = self.raw.strands, self.raw.danglings
        loop_nodes, loop_orders = self.raw.loops, self.raw.loop_orders

        refined_strands = strands.copy()
        refined_danglings = danglings.copy()

        # Step 0: crosslinkers present in the raw topology.
        crosslinkers_d3 = [node for node, degree in graph.degree if degree == 3]
        crosslinkers_d4 = [node for node, degree in graph.degree if degree == 4]

        # Step 1: crosslinkers (d3/d4) that sit on a primary loop, and
        # crosslinkers (d3) that already carry a dangling chain.
        lp_d4: list[int] = []
        lp_d3: list[int] = []
        for order, loop in zip(loop_orders, loop_nodes):
            if order != 1:
                continue
            for node in loop:
                degree = graph.degree(node)
                if degree < 3:
                    continue
                (lp_d3 if degree == 3 else lp_d4).append(node)

        dl_d3: list[int] = [
            node for dangling in danglings for node in crosslinkers_d3 if node in dangling
        ]

        for strand in strands:
            for node in lp_d3:
                if node in strand:
                    if strand in refined_strands:
                        refined_strands.remove(strand)
                    refined_danglings.append(strand)
                    if node == strand[0] and graph.degree(strand[-1]) == 3:
                        dl_d3.append(strand[-1])
                    elif node == strand[-1] and graph.degree(strand[0]) == 3:
                        dl_d3.append(strand[0])

        self._merge_strands_at(refined_strands, dl_d3)
        logger.info("Step 1 completed")

        # Step 2: crosslinkers (d4) that sit on a primary loop, or carry one
        # or two dangling chains.
        self._merge_strands_at(refined_strands, lp_d4)

        dl_d4: list[int] = []
        dl_d4_switch: list[int] = []
        for node in crosslinkers_d4:
            attached = sum(1 for dangling in danglings if node in dangling)
            if attached == 2:
                dl_d4.append(node)
            if attached == 1:
                dl_d4_switch.append(node)

        self._merge_strands_at(refined_strands, dl_d4)
        logger.info("Step 2 completed")

        # Step 3: recompute loop order, excluding crosslinkers that were
        # just reclassified as linear units (primary loops keep order 1).
        wasted_crosslinkers = lp_d4 + lp_d3 + dl_d4 + dl_d3
        refined_loop_orders: list[int] = []
        for order, loop in zip(loop_orders, loop_nodes):
            if order == 1:
                refined_loop_orders.append(order)
            else:
                refined_loop_orders.append(
                    order - sum(1 for node in loop if node in wasted_crosslinkers)
                )
        logger.info("Step 3 completed")

        # Step 4: effective vs. total crosslinker counts.
        all_crosslinker_d4 = len(crosslinkers_d4)
        all_crosslinker_d3 = len(crosslinkers_d3)
        effective_crosslinker_d4 = (
            all_crosslinker_d4 - len(lp_d4) - len(dl_d4_switch) - len(dl_d4)
        )
        effective_crosslinker_d3 = (
            all_crosslinker_d3 - len(lp_d3) - len(dl_d3) + len(dl_d4_switch)
        )
        correct_crosslinker_d4 = len(lp_d4) + len(dl_d4) + len(dl_d4_switch)
        correct_crosslinker_d3 = len(lp_d3) + len(dl_d3) - len(dl_d4_switch)
        logger.info("Step 4 completed")

        return ScanlanCaseResult(
            refined=RefinedStats(
                strands=refined_strands,
                danglings=refined_danglings,
                loops=loop_nodes,
                loop_orders=refined_loop_orders,
                effective_crosslinker_d4=effective_crosslinker_d4,
                effective_crosslinker_d3=effective_crosslinker_d3,
                all_crosslinker_d4=all_crosslinker_d4,
                all_crosslinker_d3=all_crosslinker_d3,
                correct_crosslinker_d4=correct_crosslinker_d4,
                correct_crosslinker_d3=correct_crosslinker_d3,
                switched_crosslinker_d4=len(dl_d4_switch),
            ),
            raw=self.raw,
        )

    @staticmethod
    def _merge_strands_at(strands: list[Strand], nodes: list[int]) -> None:
        """In place: splice together the (at most two) strands ending at each node."""
        for node in nodes:
            need_to_merge = []
            for strand in list(strands):
                if node == strand[0] or node == strand[-1]:
                    strands.remove(strand)
                    need_to_merge.append(strand)
            if len(need_to_merge) == 2:
                strands.append(merge_strands(node, *need_to_merge))
