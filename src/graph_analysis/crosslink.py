"""Shortest simple paths between crosslinker nodes, used to identify
elastically effective strands (crosslinker-to-crosslinker) and dangling
chains (crosslinker-to-chain-end).
"""

from __future__ import annotations

import itertools
import logging
import multiprocessing as mp
import pickle
import time
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable, Optional

import networkx as nx
import tqdm
from networkx.algorithms.simple_paths import PathBuffer, _bidirectional_shortest_path

from .network import PolymerNetwork

logger = logging.getLogger(__name__)

NodeFilter = Callable[[int], bool]

_is_crosslinker: NodeFilter = lambda degree: degree >= 3
_is_chain_end: NodeFilter = lambda degree: degree == 1


class CrosslinkAction(str, Enum):
    """Which pair of node sets to search paths between.

    CROSSLINKERS looks for strands connecting two crosslinkers (degree >= 3),
    excluding a node from being paired with itself. DANGLINGS looks for
    dangling chains connecting a crosslinker to a chain end (degree == 1),
    and keeps crosslinker neighbours in the search space since a dangling
    chain may be a single bond.
    """

    CROSSLINKERS = "crosslinkers"
    DANGLINGS = "danglings"

    @property
    def source_filter(self) -> NodeFilter:
        return _is_crosslinker

    @property
    def dest_filter(self) -> NodeFilter:
        return _is_crosslinker if self is CrosslinkAction.CROSSLINKERS else _is_chain_end

    @property
    def ignored_filter(self) -> NodeFilter:
        return _is_crosslinker

    @property
    def keep_neighbours(self) -> bool:
        return self is CrosslinkAction.DANGLINGS


@dataclass(slots=True)
class CrosslinkPathAnalyzer:
    """Finds all shortest simple paths between crosslinker nodes and a target set.

    The network is split into connected components; for each component a
    search matrix of (source, dest) node pairs is built and cached, and the
    shortest-simple-paths search for each pair is farmed out to a process
    pool.
    """

    sys_name: str
    action: CrosslinkAction
    cache_dir: Optional[Path] = None
    ncpus: Optional[int] = None
    max_chunk_size: Optional[int] = None

    def __post_init__(self) -> None:
        if self.cache_dir is None:
            self.cache_dir = Path(f"cache-{self.sys_name}")
        self.cache_dir.mkdir(exist_ok=True)

    def run(self, network: PolymerNetwork) -> list[list[int]]:
        """Run the search over every connected component of `network`."""
        components = list(network.connected_components())
        logger.info("The total number of clusters are %d", len(components))

        all_paths: list[list[int]] = []
        for idx, component in enumerate(components):
            all_paths.extend(self._run_component(idx, component))
        return all_paths

    def _run_component(self, idx: int, component: PolymerNetwork) -> list[list[int]]:
        source_nodes = component.nodes_with_degree(self.action.source_filter)
        dest_nodes = component.nodes_with_degree(self.action.dest_filter)
        ignored_nodes = component.nodes_with_degree(self.action.ignored_filter)

        if not source_nodes or not dest_nodes:
            return []

        logger.info(
            "Number of source and target nodes: %d and %d", len(source_nodes), len(dest_nodes)
        )

        prefix = f"{self.sys_name}-{idx}"
        subgraph_cache = self.cache_dir / f"{prefix}-subgraph.pkl"
        matrix_cache = self.cache_dir / f"{prefix}-crosslink_matrix.pkl"
        ignored_cache = self.cache_dir / f"{prefix}-ignored_nodes.pkl"

        with subgraph_cache.open("wb") as out:
            pickle.dump(component.graph, out)

        crosslink_matrix = self._crosslink_matrix(source_nodes, dest_nodes, matrix_cache)
        logger.info("Number of crosslinker pairs: %d", len(crosslink_matrix))

        with ignored_cache.open("wb") as out:
            pickle.dump(set(map(int, ignored_nodes)), out)

        return self._search(crosslink_matrix, subgraph_cache, ignored_cache)

    def _crosslink_matrix(
        self, source_nodes: list[int], dest_nodes: list[int], cache_file: Path
    ) -> set[tuple[int, int]]:
        """The search space: (source, dest) node pairs to find paths between.

        When the source and dest pools are the same nodes (crosslinker-to-
        crosslinker), each unordered pair is generated exactly once via
        `itertools.combinations` - equivalent to the general cross-product-
        then-dedup below, but without ever forming the full `|source| x
        |dest|` product (quadratic in the crosslinker count, which for a
        densely crosslinked network can be tens of millions of pairs before
        dedup even runs). Sorting `source_nodes` first ensures each pair
        still comes out (min, max) - matching the general branch's
        `tuple(sorted(...))` convention - since when two crosslinkers have
        more than one equally-short path between them, which one is found
        first depends on which side of the pair is passed as the search's
        "source" vs. "target".
        """
        if cache_file.exists():
            with cache_file.open("rb") as inp:
                return pickle.load(inp)

        keep_neighbours = self.action.keep_neighbours
        if not keep_neighbours and source_nodes == dest_nodes:
            matrix = set(itertools.combinations(sorted(map(int, source_nodes)), 2))
        else:
            matrix = {
                tuple(sorted((int(s), int(d))))
                for s in source_nodes
                for d in dest_nodes
                if keep_neighbours or s != d
            }
        with cache_file.open("wb") as out:
            pickle.dump(matrix, out)
        return matrix

    def _search(
        self, crosslink_matrix: set[tuple[int, int]], subgraph_cache: Path, ignored_cache: Path
    ) -> list[list[int]]:
        num_workers = max(self.ncpus or (mp.cpu_count() - 2), 1)
        chunk_size = max(len(crosslink_matrix) // num_workers, 1)
        if self.max_chunk_size:
            chunk_size = min(chunk_size, self.max_chunk_size)

        logger.info("Number of workers: %d, chunk_size: %d", num_workers, chunk_size)

        start = time.time()
        # Extend a single flat list as results arrive, rather than collecting
        # one (mostly-empty) list per pair and flattening at the end: for a
        # densely crosslinked network the pair count can be tens of millions,
        # and holding a separate result object per pair (plus a second,
        # flattened copy) is the dominant memory cost.
        all_paths: list[list[int]] = []
        with mp.Pool(num_workers, _init_worker, (subgraph_cache, ignored_cache)) as pool:
            for paths in tqdm.tqdm(
                pool.imap(_all_shortest_paths, crosslink_matrix, chunksize=chunk_size),
                total=len(crosslink_matrix),
            ):
                all_paths.extend(paths)
        logger.info("Total time %.2fs", time.time() - start)
        return all_paths


# --- multiprocessing worker internals ---------------------------------------
#
# These stay as plain module-level functions (rather than methods) so they
# pickle cleanly under the "spawn" start method that macOS and Windows use by
# default: each worker process re-imports this module, calls `_init_worker`
# once to populate its own copy of these globals from the cache files, and
# then repeatedly calls `_all_shortest_paths`.

_worker_graph: Optional[nx.Graph] = None
_worker_ignored_nodes: Optional[set[int]] = None


def _init_worker(subgraph_cache: Path, ignored_cache: Path) -> None:
    global _worker_graph, _worker_ignored_nodes
    with open(subgraph_cache, "rb") as inp:
        _worker_graph = pickle.load(inp)
    with open(ignored_cache, "rb") as inp:
        _worker_ignored_nodes = pickle.load(inp)


def _shortest_simple_paths(graph: nx.Graph, source: int, target: int, ignored: set[int]):
    """Generate simple paths from source to target, shortest first.

    Adapted from `networkx.shortest_simple_paths` to support ignoring a set
    of nodes (crosslinkers other than source/target) throughout the search -
    a parameter the public networkx API doesn't expose.
    """
    if source not in graph:
        raise nx.NodeNotFound(f"source node {source} not in graph")
    if target not in graph:
        raise nx.NodeNotFound(f"target node {target} not in graph")

    ignored = ignored - {source, target}

    list_a: list[list[int]] = []
    list_b = PathBuffer()
    prev_path = None
    while True:
        if not prev_path:
            length, path = _bidirectional_shortest_path(
                graph, source, target, ignore_nodes=ignored
            )
            list_b.push(length, path)
        else:
            ignore_nodes = ignored
            ignore_edges = set()
            for i in range(1, len(prev_path)):
                root = prev_path[:i]
                root_length = len(root)
                for path in list_a:
                    if path[:i] == root:
                        ignore_edges.add((path[i - 1], path[i]))
                try:
                    length, spur = _bidirectional_shortest_path(
                        graph,
                        root[-1],
                        target,
                        ignore_nodes=ignore_nodes,
                        ignore_edges=ignore_edges,
                    )
                    path = root[:-1] + spur
                    list_b.push(root_length + length, path)
                except nx.NetworkXNoPath:
                    pass
                ignore_nodes.add(root[-1])

        if list_b:
            path = list_b.pop()
            yield path
            list_a.append(path)
            prev_path = path
        else:
            break


def _all_shortest_paths(source_target: tuple[int, int]) -> list[list[int]]:
    """Worker task: every shortest simple path between one (source, target) pair."""
    source, target = source_target
    paths: list[list[int]] = []
    try:
        for path in _shortest_simple_paths(_worker_graph, source, target, _worker_ignored_nodes):
            if path:
                paths.append(path)
    except nx.NetworkXNoPath:
        pass
    return paths
