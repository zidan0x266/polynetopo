"""End-to-end network analysis in a single pass.

Loads a topology once and runs every downstream stage against that one
in-memory graph - strands, dangling chains, minimum-cycle-basis loops, and
their Scanlan-Case refinement - instead of writing/reading an intermediate
pickle between each step.
"""

from __future__ import annotations

import logging
from typing import Optional

from .crosslink import CrosslinkAction, CrosslinkPathAnalyzer
from .loops import LoopAnalyzer
from .network import PolymerNetwork
from .scanlan_case import NetworkStats, ScanlanCaseAnalyzer, ScanlanCaseResult

logger = logging.getLogger(__name__)


class NetworkAnalysisPipeline:
    """Runs strands + danglings + loops + Scanlan-Case from one loaded network."""

    def __init__(self, network: PolymerNetwork, n_vertices: Optional[int] = None) -> None:
        self.network = network
        self.n_vertices = n_vertices

    def run(self) -> ScanlanCaseResult:
        raw = self.raw_stats()
        return ScanlanCaseAnalyzer(self.network.graph, raw).run()

    def raw_stats(self) -> NetworkStats:
        strands = self._search(CrosslinkAction.CROSSLINKERS)
        danglings = self._search(CrosslinkAction.DANGLINGS)

        loops = LoopAnalyzer.from_polymer_network(self.network, self.n_vertices).minimum_cycle_basis()
        logger.info("Found %d strands, %d danglings, %d loops", len(strands), len(danglings), len(loops))

        return NetworkStats.from_loop_basis(strands, danglings, loops)

    def _search(self, action: CrosslinkAction) -> list[list[int]]:
        return CrosslinkPathAnalyzer(action).run(self.network)
