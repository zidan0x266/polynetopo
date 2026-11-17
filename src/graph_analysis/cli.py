"""Command-line entry point for the polymer network analysis pipeline.

A single `--function` flag selects what to run:

    topology       strands, danglings, loops and their Scanlan-Case refinement
                    (--stage picks full/strands/danglings/loops/scanlan-case; default: full)
    shortest-path  shortest paths across the periodic boundary (Yu & Jackson 2025),
                    which govern where and when the first bond breaks
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
from .pipeline import NetworkAnalysisPipeline
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


def _sibling_path(output_file: str, suffix: str) -> Path:
    """`results.pkl` -> `results_<suffix>.pkl`."""
    path = Path(output_file)
    return path.with_stem(f"{path.stem}_{suffix}")


def _dump(obj: object, path: Path | str) -> None:
    with open(path, "wb") as out:
        pickle.dump(obj, out)


def _write_result(
    result: ScanlanCaseResult, output_file: str, write_raw: bool, write_stage_files: bool = False
) -> None:
    _dump(result.refined, output_file)
    logger.info("Wrote refined result to %s", output_file)

    if write_stage_files:
        # The inputs of `--stage scanlan-case`, so it can be rerun without the slow loop search.
        raw = result.raw
        for suffix, content in (
            ("strands", raw.strands),
            ("danglings", raw.danglings),
            ("loops", [raw.loops, raw.loop_orders]),
        ):
            stage_file = _sibling_path(output_file, suffix)
            _dump(content, stage_file)
            logger.info("Wrote %s to %s", suffix, stage_file)

    if write_raw:
        raw_file = _sibling_path(output_file, "raw")
        _dump(result.raw, raw_file)
        logger.info("Wrote raw (pre-refinement) result to %s", raw_file)


# --- topology ----------------------------------------------------------------


def _topology_full(args: argparse.Namespace) -> None:
    network = _load_network(args)

    pipeline = NetworkAnalysisPipeline(network=network, n_vertices=args.nvertices)
    result = pipeline.run()
    _write_result(result, args.output_file, args.raw, write_stage_files=True)


def _topology_paths(action: CrosslinkAction):
    """A stage that writes just one flat list of strand or dangling paths."""

    def stage(args: argparse.Namespace) -> None:
        network = _load_network(args)
        paths = CrosslinkPathAnalyzer(action).run(network)
        _dump(paths, args.output_file)
        logger.info("Wrote %d %s to %s", len(paths), action.value, args.output_file)

    return stage


_topology_strands = _topology_paths(CrosslinkAction.CROSSLINKERS)
_topology_danglings = _topology_paths(CrosslinkAction.DANGLINGS)


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
    "strands": _topology_strands,
    "danglings": _topology_danglings,
    "loops": _topology_loops,
    "scanlan-case": _topology_scanlan_case,
}


def _topology_outputs(args: argparse.Namespace) -> list[Path]:
    """Every file a topology stage writes."""
    outputs = [Path(args.output_file)]
    if args.stage == "full":
        outputs += [_sibling_path(args.output_file, suffix) for suffix in ("strands", "danglings", "loops")]
    if args.raw and args.stage in ("full", "scanlan-case"):
        outputs.append(_sibling_path(args.output_file, "raw"))
    return outputs


def _topology(args: argparse.Namespace) -> None:
    existing = [str(path) for path in _topology_outputs(args) if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(f"Output already exists; use --overwrite or --output-file: {', '.join(existing)}")
    _TOPOLOGY_STAGES[args.stage](args)


# --- shortest-path -------------------------------------------------------------


def _shortest_path(args: argparse.Namespace) -> None:
    from .lammps import Direction
    from .shortest_path import PeriodicShortestPathAnalysis, ShortestPathConfig

    config = ShortestPathConfig(
        direction=Direction.parse(args.direction),
        max_period=args.max_period,
        buffer_images=args.buffer_images,
        slice_max=args.slice_max,
        remove_redundant=not args.keep_redundant,
        min_sp_only=args.min_sp_only,
        ncpus=args.ncpus,
        checkpoint_file=args.checkpoint_file,
    )
    analysis = PeriodicShortestPathAnalysis.from_lammps(args.input_file, config)

    suffix = f".{config.direction.name.lower()}.sp.pkl"
    output = (
        Path(args.output_file) if args.output_file else Path(args.input_file).with_suffix(suffix)
    ).expanduser().resolve()
    if not args.overwrite and output.exists():
        raise FileExistsError(f"Output already exists; use --overwrite or --output-file: {output}")

    result = analysis.run()
    archive = result.save(output, overwrite=args.overwrite)
    print(json.dumps({**result.summary(), "output": str(archive)}, indent=2, allow_nan=False))


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

    from .betweenness import AnalysisConfig, SpatialBetweennessAnalysis
    from .lammps import Direction

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

    suffix = f".{config.direction.name.lower()}.spatial.pkl"
    if args.sample_vertices is not None:
        suffix = f".sample{args.sample_vertices}" + suffix
    output = (
        Path(args.output_file) if args.output_file else Path(args.input_file).with_suffix(suffix)
    ).expanduser().resolve()
    if output.suffix != ".pkl":
        raise ValueError("Output filename must end in .pkl")
    if analysis.network.source is not None and Path(analysis.network.source) == output:
        raise ValueError("Output file must not replace the input LAMMPS file")
    if not args.overwrite and output.exists():
        raise FileExistsError(f"Output already exists; use --overwrite or --output-file: {output}")

    if not args.quiet:
        ig.set_progress_handler(_ProgressLogger())
    try:
        result = analysis.run()
    finally:
        if not args.quiet:
            ig.set_progress_handler(None)

    archive = result.save(output, overwrite=args.overwrite)
    print(json.dumps({**result.summary(), "output": str(archive)}, indent=2, allow_nan=False))


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
    ("topology", "strands"): ["input_file", "output_file"],
    ("topology", "danglings"): ["input_file", "output_file"],
    ("topology", "loops"): ["input_file", "output_file"],
    ("topology", "scanlan-case"): [
        "input_file", "strands_file", "danglings_file", "loops_file", "output_file",
    ],
    ("shortest-path", None): ["input_file"],
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
        "polynetopo",
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
        help="Topology file format for --function topology: a full LAMMPS "
        "data file (default), a pre-trimmed bare bond list, or an RMC-style edge list",
    )
    common.add_argument(
        "--nvertices", type=int, default=None,
        help="(topology stage=full/loops) Largest atom ID (the loop search uses atom IDs as "
        "vertex IDs, so this is not the number of atoms); inferred from the topology if omitted",
    )

    topology = parser.add_argument_group("--function topology")
    topology.add_argument(
        "--stage",
        choices=list(_TOPOLOGY_STAGES),
        default="full",
        help="full (default): strands, danglings, loops and Scanlan-Case in one pass, also "
        "writing <output>_strands/_danglings/_loops.pkl. strands/danglings/loops/scanlan-case: "
        "just that stage (scanlan-case needs the three files from the others), e.g. to rerun "
        "the refinement or run the memory-heavy loop search on a bigger machine separately",
    )
    topology.add_argument(
        "--raw",
        action="store_true",
        help="(stage full/scanlan-case) Also write the pre-refinement strands/danglings/loops "
        "to <output>_raw.<ext>",
    )
    topology.add_argument("--strands-file", help="(stage scanlan-case) Pickled strands, from --stage strands or <output>_strands.pkl")
    topology.add_argument("--danglings-file", help="(stage scanlan-case) Pickled danglings, from --stage danglings or <output>_danglings.pkl")
    topology.add_argument("--loops-file", help="(stage scanlan-case) Pickled loop basis, from --stage loops or <output>_loops.pkl")

    common.add_argument(
        "--direction", default="x", metavar="{x,y,z,0,1,2}",
        help="(betweenness) loading axis; (shortest-path) periodic axis the paths wrap along (default: x)",
    )
    common.add_argument(
        "--overwrite", action="store_true", help="Replace existing output files (default: stop if one exists)"
    )

    shortest_path = parser.add_argument_group("--function shortest-path")
    shortest_path.add_argument(
        "--max-period", type=int, default=8, metavar="M",
        help="The paper's M: search periods p = 1..M, each wrapping p box lengths (default: 8)",
    )
    shortest_path.add_argument(
        "--buffer-images", type=int, default=1,
        help="Extra periodic images on each side, so paths can leave the stack between their "
        "ends (default: 1; 0 keeps them inside it)",
    )
    shortest_path.add_argument(
        "--slice-max", type=float, metavar="L",
        help="Only start paths from atoms within L of the box's low face along --direction "
        "(default: every atom); reduces cost, and keeps only the SPs found from those atoms",
    )
    shortest_path.add_argument(
        "--keep-redundant", action="store_true", help="Skip the removal of redundant paths (Cases 1 and 2)"
    )
    shortest_path.add_argument(
        "--min-sp-only", action="store_true",
        help="Also apply the paper's optional Case 3: keep only the shortest SP of each start atom",
    )
    shortest_path.add_argument(
        "--ncpus", type=int, default=None, help="Worker processes for the per-atom search (default: 1)"
    )
    shortest_path.add_argument(
        "--checkpoint-file", type=Path,
        help="Save progress after each period here and resume from it if it exists; "
        "deleted when the search completes",
    )

    betweenness = parser.add_argument_group("--function betweenness")
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
