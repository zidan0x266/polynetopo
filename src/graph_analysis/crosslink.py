"""Chains between crosslinker nodes, used to identify elastically effective
strands (crosslinker-to-crosslinker) and dangling chains (crosslinker-to-
chain-end).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Callable

from .network import PolymerNetwork

logger = logging.getLogger(__name__)

NodeFilter = Callable[[int], bool]

_is_crosslinker: NodeFilter = lambda degree: degree >= 3
_is_chain_end: NodeFilter = lambda degree: degree == 1


class CrosslinkAction(str, Enum):
    """Which pair of node sets to find chains between.

    CROSSLINKERS looks for strands connecting two different crosslinkers
    (degree >= 3). DANGLINGS looks for dangling chains connecting a
    crosslinker to a chain end (degree == 1); a dangling chain may be a
    single bond.
    """

    CROSSLINKERS = "crosslinkers"
    DANGLINGS = "danglings"

    @property
    def source_filter(self) -> NodeFilter:
        return _is_crosslinker

    @property
    def dest_filter(self) -> NodeFilter:
        return _is_crosslinker if self is CrosslinkAction.CROSSLINKERS else _is_chain_end


@dataclass(slots=True)
class CrosslinkPathAnalyzer:
    """Finds every chain between crosslinkers and a target set.

    Every path of interest is a chain of degree-2 nodes (every other
    crosslinker is ignored) from a crosslinker to either another crosslinker
    (strand) or a degree-1 chain end (dangling). So instead of searching each
    (source, dest) pair, this follows each bond leaving a crosslinker until the
    chain hits a node that isn't degree 2: one O(edges) pass with no pair
    matrix, process pool or per-worker graph copies. Paths are oriented from
    the lower to the higher node ID.

    This finds every distinct chain between a crosslinker pair, including when
    three or more parallel chains join the same two crosslinkers. (Releases
    before the chain walk used a pairwise Yen's-algorithm search that missed
    some of those, about 0.3% of strands on the 104,704-atom sys3 network.)
    """

    action: CrosslinkAction

    def run(self, network: PolymerNetwork) -> list[list[int]]:
        graph = network.graph
        adj = graph.adj
        degree = dict(graph.degree())
        is_source = self.action.source_filter
        is_dest = self.action.dest_filter
        strands = self.action is CrosslinkAction.CROSSLINKERS

        paths: list[list[int]] = []
        for source in sorted(n for n, d in degree.items() if is_source(d)):
            for first in sorted(adj[source]):
                path = [source, first]
                prev, node = source, first
                while degree[node] == 2:
                    nxt = next(n for n in adj[node] if n != prev)
                    prev, node = node, nxt
                    path.append(node)
                if not is_dest(degree[node]) or (strands and node == source):
                    continue
                if strands and node < source:
                    continue  # the same strand is found from its other end
                if node < source:
                    path.reverse()
                paths.append(path)
        logger.info("Found %d %s paths", len(paths), self.action.value)
        return paths
