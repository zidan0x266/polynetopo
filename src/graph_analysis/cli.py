"""Command-line entry point for the polymer network analysis pipeline.

A single `--function` flag selects what to run:

    topology       strands, danglings, loops and their Scanlan-Case refinement
                    (--stage picks full/loops/scanlan-case; default: full)
    shortest-path  shortest simple paths between crosslinkers, or crosslinker
                    to chain-end (dangling chains)
    betweenness    native spatially-weighted edge betweenness (GEBC/OBC/DBC)
    compare        print the entry counts of two pickled result files
"""

from __future__ import annotations

import argparse
import json
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


_FORMAT_LOADERS = {
    "lammps": PolymerNetwork.from_lammps,
    "bond-list": PolymerNetwork.from_bond_list,
    "rmc": PolymerNetwork.from_rmc,
}


def _load_network(args: argparse.Namespace) -> PolymerNetwork:
    network = _FORMAT_LOADERS[args.format](args.input_file)
    logger.info("Loaded %s", network)
    return network


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


# --- topology ----------------------------------------------------------------


def _topology_full(args: argparse.Namespace) -> None:
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


def _topology_loops(args: argparse.Namespace) -> None:
    network = _load_network(args)
    analyzer = LoopAnalyzer.from_polymer_network(network, args.nvertices)
    logger.info(
        "Graph has %d nodes and %d edges", analyzer.graph.vcount(), analyzer.graph.ecount()
    )
    basis = analyzer.minimum_cycle_basis()

    with open(args.output_file, "wb") as out:
        pickle.dump(basis.as_legacy_tuple(), out)
    logger.info("Wrote %d loops to %s", len(basis), args.output_file)


def _topology_scanlan_case(args: argparse.Namespace) -> None:
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


_TOPOLOGY_STAGES = {
    "full": _topology_full,
    "loops": _topology_loops,
    "scanlan-case": _topology_scanlan_case,
}


def _topology(args: argparse.Namespace) -> None:
    _TOPOLOGY_STAGES[args.stage](args)


# --- shortest-path -------------------------------------------------------------


def _shortest_path(args: argparse.Namespace) -> None:
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


# --- betweenness ---------------------------------------------------------------


class _ProgressLogger:
    """Throttle igraph's per-source progress messages."""

    def __init__(self) -> None:
        self._message = ""
        self._percent = -10.0

    def __call__(self, message: str, percent: float) -> None:
        if message != self._message or percent >= self._percent + 5 or percent == 100:
            logger.info("%s%.0f%%", message, percent)
            self._message, self._percent = message, percent


def _betweenness(args: argparse.Namespace) -> None:
    # Imported lazily: this subpackage needs the in-house igraph build (see
    # vendor/README.md) only at analysis time, not merely to parse arguments.
    import igraph as ig

    from .betweenness import AnalysisConfig, Direction, SpatialBetweennessAnalysis

    if not args.quiet:
        logging.getLogger().setLevel(logging.INFO)

    config = AnalysisConfig(
        direction=Direction.parse(args.direction),
        validate_stock=args.validate_stock,
        dbc_length=args.dbc_length,
    )
    analysis = SpatialBetweennessAnalysis.from_lammps(args.input_file, config)

    if args.sample_vertices is not None:
        analysis = analysis.sample(args.sample_vertices)
        logger.warning("Analyzing an induced debug subgraph; results are not full-network centralities")

    if args.inspect:
        print(json.dumps(analysis.inspect(), indent=2, allow_nan=False))
        return

    suffix = f".{config.direction.name.lower()}.spatial.npz"
    if args.sample_vertices is not None:
        suffix = f".sample{args.sample_vertices}" + suffix
    output = (
        Path(args.output_file) if args.output_file else Path(args.input_file).with_suffix(suffix)
    ).expanduser().resolve()
    if output.suffix != ".npz":
        raise ValueError("Output filename must end in .npz")
    if analysis.network.source in (output, output.with_suffix(".json")):
        raise ValueError("Output files must not replace the input LAMMPS file")
    if not args.overwrite and (output.exists() or output.with_suffix(".json").exists()):
        raise FileExistsError(f"Output already exists; use --overwrite or --output-file: {output}")

    if not args.quiet:
        ig.set_progress_handler(_ProgressLogger())
    try:
        result = analysis.run()
    finally:
        if not args.quiet:
            ig.set_progress_handler(None)

    archive, report = result.save(output, overwrite=args.overwrite)
    print(json.dumps({**result.summary(), "output": str(archive), "report": str(report)}, indent=2, allow_nan=False))


# --- compare -------------------------------------------------------------------


def _compare(args: argparse.Namespace) -> None:
    with open(args.file1, "rb") as inp:
        data1 = pickle.load(inp)
    with open(args.file2, "rb") as inp:
        data2 = pickle.load(inp)
    print(len(data1), len(data2))


# --- argument parsing and dispatch ----------------------------------------------

_DISPATCH = {
    "topology": _topology,
    "shortest-path": _shortest_path,
    "betweenness": _betweenness,
    "compare": _compare,
}

# (function, stage) -> required argument names; stage is only meaningful for
# "topology" and is otherwise ignored.
_REQUIRED_ARGS = {
    ("topology", "full"): ["input_file", "output_file"],
    ("topology", "loops"): ["input_file", "output_file"],
    ("topology", "scanlan-case"): [
        "input_file", "strands_file", "danglings_file", "loops_file", "output_file",
    ],
    ("shortest-path", None): ["sys_name", "input_file", "output_file", "action"],
    ("betweenness", None): ["input_file"],
    ("compare", None): ["file1", "file2"],
}


def _validate(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    stage = args.stage if args.function == "topology" else None
    required = _REQUIRED_ARGS[(args.function, stage)]
    missing = [name for name in required if getattr(args, name, None) is None]
    if missing:
        flags = ", ".join(f"--{name.replace('_', '-')}" for name in missing)
        where = f"--function {args.function}" + (f" --stage {stage}" if stage else "")
        parser.error(f"{where} requires: {flags}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        "graph-analysis",
        description="Topology, shortest-path and spatial betweenness analysis for crosslinked polymer networks",
    )
    parser.add_argument(
        "--function",
        choices=list(_DISPATCH),
        default="topology",
        help="Which analysis to run (default: topology)",
    )

    common = parser.add_argument_group("shared")
    common.add_argument("--input-file", help="Topology file (topology/shortest-path/betweenness)")
    common.add_argument("--output-file", help="Output file; see each function's section below for the default")
    common.add_argument(
        "--format",
        choices=list(_FORMAT_LOADERS),
        default="lammps",
        help="Topology file format for --function topology/shortest-path: a full LAMMPS "
        "data file (default), a pre-trimmed bare bond list, or an RMC-style edge list",
    )
    common.add_argument(
        "--sys-name", default=None,
        help="(topology stage=full, shortest-path) Cache directory prefix; for topology, "
        "defaults to the output file's stem",
    )
    common.add_argument(
        "--nvertices", type=int, default=None,
        help="(topology stage=full/loops) Number of nodes in the system; inferred from the "
        "topology if omitted",
    )
    common.add_argument("--ncpus", type=int, default=None, help="(topology stage=full, shortest-path)")
    common.add_argument("--max-chunk-size", type=int, default=None, help="(topology stage=full, shortest-path)")

    topology = parser.add_argument_group("--function topology")
    topology.add_argument(
        "--stage",
        choices=list(_TOPOLOGY_STAGES),
        default="full",
        help="full: strands+danglings+loops+Scanlan-Case in one pass (default). "
        "loops/scanlan-case: just that stage, e.g. to resume after a crash or run "
        "a memory-heavy stage on a bigger machine separately",
    )
    topology.add_argument(
        "--raw",
        action="store_true",
        help="(stage full/scanlan-case) Also write the pre-refinement strands/danglings/loops "
        "to <output>_raw.<ext>",
    )
    topology.add_argument("--strands-file", help="(stage scanlan-case) Pickled strands, from --function shortest-path --action crosslinkers")
    topology.add_argument("--danglings-file", help="(stage scanlan-case) Pickled danglings, from --function shortest-path --action danglings")
    topology.add_argument("--loops-file", help="(stage scanlan-case) Pickled loop basis, from --stage loops")

    shortest_path = parser.add_argument_group("--function shortest-path")
    shortest_path.add_argument("--action", choices=[a.value for a in CrosslinkAction])

    betweenness = parser.add_argument_group("--function betweenness")
    betweenness.add_argument(
        "--direction", default="x", metavar="{x,y,z,0,1,2}", help="Loading axis (default: x)"
    )
    betweenness.add_argument(
        "--dbc-length", type=float, metavar="L",
        help="Optional length divisor for postprocessed DBC only; raw DBC always uses distance alone",
    )
    betweenness.add_argument(
        "--validate-stock", action="store_true", help="Run a separate stock GEBC pass and require agreement"
    )
    betweenness.add_argument(
        "--inspect", action="store_true", help="Validate input and print topology without computing centrality"
    )
    betweenness.add_argument(
        "--sample-vertices", type=int, metavar="N",
        help="Debug only: induced neighborhood, not full-network centrality",
    )
    betweenness.add_argument(
        "--overwrite", action="store_true", help="Replace existing output files"
    )
    betweenness.add_argument("--quiet", action="store_true", help="Suppress progress logging")

    compare = parser.add_argument_group("--function compare")
    compare.add_argument("--file1")
    compare.add_argument("--file2")

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    _validate(parser, args)
    _DISPATCH[args.function](args)


if __name__ == "__main__":
    main()
