"""Command-line entry point for the polymer network analysis pipeline.

Subcommands:
    analyze       strands, danglings, loops and their Scanlan-Case refinement,
                  all in one pass over a single topology load
    distances     shortest simple paths between crosslinkers, or crosslinker
                  to chain-end (dangling chains) - one stage of `analyze`
    loops         minimum cycle basis loop detection - one stage of `analyze`
    scanlan-case  Scanlan-Case crosslinker density refinement, from files
                  produced by `distances`/`loops`
    compare       print the entry counts of two pickled result files
"""

from __future__ import annotations

import argparse
import logging
import pickle
from pathlib import Path

from .crosslink import CrosslinkAction, CrosslinkPathAnalyzer
from .loops import LoopAnalyzer, LoopBasis
from .network import PolymerNetwork
from .pipeline import NetworkAnalysisPipeline, PipelineConfig
from .scanlan_case import NetworkStats, ScanlanCaseAnalyzer, ScanlanCaseResult

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("graph_analysis")


_TOPOLOGY_LOADERS = {
    "lammps": PolymerNetwork.from_lammps,
    "bond-list": PolymerNetwork.from_bond_list,
    "rmc": PolymerNetwork.from_rmc,
}


def _load_network(args: argparse.Namespace) -> PolymerNetwork:
    network = _TOPOLOGY_LOADERS[args.topology](args.input_file)
    logger.info("Loaded %s", network)
    return network


def _add_topology_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--topology",
        choices=list(_TOPOLOGY_LOADERS),
        default="lammps",
        help="Topology file format: a full LAMMPS data file (default), a "
        "pre-trimmed bare bond list, or an RMC-style edge list",
    )


def _raw_output_path(output_file: str) -> Path:
    """`results.pkl` -> `results_raw.pkl`."""
    path = Path(output_file)
    return path.with_stem(path.stem + "_raw")


def _write_result(result: ScanlanCaseResult, output_file: str, write_raw: bool) -> None:
    with open(output_file, "wb") as out:
        pickle.dump(result.refined, out)
    logger.info("Wrote refined result to %s", output_file)

    if write_raw:
        raw_file = _raw_output_path(output_file)
        with open(raw_file, "wb") as out:
            pickle.dump(result.raw, out)
        logger.info("Wrote raw (pre-refinement) result to %s", raw_file)


def _analyze(args: argparse.Namespace) -> None:
    network = _load_network(args)

    pipeline = NetworkAnalysisPipeline(
        network=network,
        config=PipelineConfig(
            sys_name=args.sys_name or Path(args.output_file).stem,
            ncpus=args.ncpus,
            max_chunk_size=args.max_chunk_size,
        ),
        n_vertices=args.nvertices,
    )
    result = pipeline.run()
    _write_result(result, args.output_file, args.raw)


def _distances(args: argparse.Namespace) -> None:
    network = _load_network(args)

    analyzer = CrosslinkPathAnalyzer(
        sys_name=args.sys_name,
        action=CrosslinkAction(args.action),
        ncpus=args.ncpus,
        max_chunk_size=args.max_chunk_size,
    )
    paths = [path for path in analyzer.run(network) if path]

    with open(args.output_file, "wb") as out:
        pickle.dump(paths, out)
    logger.info("Wrote %d paths to %s", len(paths), args.output_file)


def _loops(args: argparse.Namespace) -> None:
    network = _load_network(args)
    analyzer = LoopAnalyzer.from_polymer_network(network, args.nvertices)
    logger.info(
        "Graph has %d nodes and %d edges", analyzer.graph.vcount(), analyzer.graph.ecount()
    )
    basis = analyzer.minimum_cycle_basis()

    with open(args.output_file, "wb") as out:
        pickle.dump(basis.as_legacy_tuple(), out)
    logger.info("Wrote %d loops to %s", len(basis), args.output_file)


def _scanlan_case(args: argparse.Namespace) -> None:
    network = _load_network(args)

    with open(args.strands_file, "rb") as inp:
        strands = pickle.load(inp)
    with open(args.danglings_file, "rb") as inp:
        danglings = pickle.load(inp)
    with open(args.loops_file, "rb") as inp:
        loops = LoopBasis(*pickle.load(inp))

    raw = NetworkStats.from_loop_basis(strands, danglings, loops)
    result = ScanlanCaseAnalyzer(network.graph, raw).run()
    _write_result(result, args.output_file, args.raw)


def _compare(args: argparse.Namespace) -> None:
    with open(args.file1, "rb") as inp:
        data1 = pickle.load(inp)
    with open(args.file2, "rb") as inp:
        data2 = pickle.load(inp)
    print(len(data1), len(data2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        "graph-analysis", description="Topology analysis pipeline for crosslinked polymer networks"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser(
        "analyze",
        help="Strands, danglings, loops and Scanlan-Case refinement in one pass",
    )
    analyze.add_argument("--input-file", required=True, help="Topology edge-list file")
    analyze.add_argument("--output-file", required=True, help="Where to write the refined result")
    analyze.add_argument(
        "--raw",
        action="store_true",
        help="Also write the pre-refinement strands/danglings/loops to <output>_raw.<ext>",
    )
    analyze.add_argument(
        "--sys-name", default=None, help="Cache directory prefix; defaults to the output file's stem"
    )
    analyze.add_argument(
        "--nvertices",
        type=int,
        default=None,
        help="Number of nodes in the system; inferred from the topology if omitted",
    )
    analyze.add_argument("--ncpus", type=int, default=None)
    analyze.add_argument("--max-chunk-size", type=int, default=None)
    _add_topology_argument(analyze)
    analyze.set_defaults(func=_analyze)

    distances = subparsers.add_parser(
        "distances", help="Shortest paths between crosslinkers, or crosslinker to chain-end"
    )
    distances.add_argument("--sys-name", required=True, help="Run name; used for the cache directory")
    distances.add_argument("--input-file", required=True, help="Topology edge-list file")
    distances.add_argument("--output-file", required=True, help="Where to write the pickled path list")
    distances.add_argument("--ncpus", type=int, default=None)
    distances.add_argument("--max-chunk-size", type=int, default=None)
    _add_topology_argument(distances)
    distances.add_argument("--action", choices=[a.value for a in CrosslinkAction], required=True)
    distances.set_defaults(func=_distances)

    loops = subparsers.add_parser("loops", help="Minimum cycle basis loop detection")
    loops.add_argument("--input-file", required=True)
    loops.add_argument("--output-file", required=True)
    loops.add_argument(
        "--nvertices",
        type=int,
        default=None,
        help="Number of nodes in the system; inferred from the topology if omitted",
    )
    _add_topology_argument(loops)
    loops.set_defaults(func=_loops)

    scanlan_case = subparsers.add_parser(
        "scanlan-case", help="Scanlan-Case crosslinker density refinement"
    )
    scanlan_case.add_argument("--input-file", required=True, help="Topology edge-list file")
    scanlan_case.add_argument("--strands-file", required=True, help="Pickled strands, from `distances --action crosslinkers`")
    scanlan_case.add_argument("--danglings-file", required=True, help="Pickled danglings, from `distances --action danglings`")
    scanlan_case.add_argument("--loops-file", required=True, help="Pickled loop basis, from `loops`")
    scanlan_case.add_argument("--output-file", required=True)
    _add_topology_argument(scanlan_case)
    scanlan_case.add_argument(
        "--raw",
        action="store_true",
        help="Also write the pre-refinement strands/danglings/loops to <output>_raw.<ext>",
    )
    scanlan_case.set_defaults(func=_scanlan_case)

    compare = subparsers.add_parser("compare", help="Print the entry counts of two pickled result files")
    compare.add_argument("file1")
    compare.add_argument("file2")
    compare.set_defaults(func=_compare)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
