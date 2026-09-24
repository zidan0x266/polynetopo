#!/usr/bin/env python3
"""Compare new-code outputs for examples/shortest_path against the old ones.

Usage (after running `loops` + `scanlan-case`, or `analyze`, with the new code):

    python3 compare_new_vs_old.py --loops igloops_new.pkl --scanlan-case scanlan_case_new.pkl
    python3 compare_new_vs_old.py --analyze-output results_new.pkl --analyze-raw results_new_raw.pkl

Prints exact-match booleans for strands, danglings, loops, loop orders, and the
Scanlan-Case crosslinker-density stats. A small number of differing strands (the
same ~0.01-0.3% seen on examples/atrp) is expected and fine: it's the intentional
fix to a list-mutation-during-iteration bug in the original strand-merging code
(see src/graph_analysis/scanlan_case.py's module docstring) - everything else
should match exactly.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

HERE = Path(__file__).parent


def norm_paths(paths):
    return sorted(tuple(p) if p[0] < p[-1] else tuple(p[::-1]) for p in paths)


def load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def compare_strands_danglings_loops(strands_new, danglings_new, loops_new_tuple):
    strands_old = load(HERE / "strands.pkl")
    danglings_old = load(HERE / "danglings.pkl")
    loops_old, loop_orders_old = load(HERE / "igloops.pkl")

    print("--- raw strands/danglings/loops ---")
    print("strands   :", len(strands_new), "vs", len(strands_old),
          "-> exact match:", norm_paths(strands_new) == norm_paths(strands_old))
    print("danglings :", len(danglings_new), "vs", len(danglings_old),
          "-> exact match:", norm_paths(danglings_new) == norm_paths(danglings_old))

    loops_new, loop_orders_new = loops_new_tuple
    print("loops     :", len(loops_new), "vs", len(loops_old),
          "-> exact match:", loops_new == loops_old)
    print("loop_orders exact match:", loop_orders_new == loop_orders_old)


def compare_scanlan_case(refined, raw):
    with open(HERE / "Scanlan_Case.pkl", "rb") as f:
        old = pickle.load(f)  # [[refined...], [raw...], [stats...]]

    stats_new = [
        refined.effective_crosslinker_d4,
        refined.effective_crosslinker_d3,
        refined.all_crosslinker_d4,
        refined.all_crosslinker_d3,
        refined.correct_crosslinker_d4,
        refined.correct_crosslinker_d3,
        refined.switched_crosslinker_d4,
    ]
    print("\n--- Scanlan-Case ---")
    print("stats new:", stats_new)
    print("stats old:", old[2])
    print("stats exact match:", stats_new == old[2])

    print("refined danglings exact match:", refined.danglings == old[0][1])
    print("refined loop_orders exact match:", refined.loop_orders == old[0][3])

    a, b = norm_paths(refined.strands), norm_paths(old[0][0])
    same_len = len(a) == len(b)
    differ = sum(1 for x, y in zip(a, b) if x != y) if same_len else "n/a"
    print(f"refined strands: same count={same_len}, differing={differ} "
          f"(a small number here is expected - see module docstring)")

    if raw is not None:
        print("\n--- raw (pre-refinement), from the --raw output ---")
        print("raw strands exact match:", norm_paths(raw.strands) == norm_paths(old[1][0]))
        print("raw danglings exact match:", norm_paths(raw.danglings) == norm_paths(old[1][1]))
        print("raw loops exact match:", raw.loops == old[1][2])
        print("raw loop_orders exact match:", raw.loop_orders == old[1][3])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strands", help="new strands.pkl (defaults to the already-validated local copy)")
    parser.add_argument("--danglings", help="new danglings.pkl (defaults to the already-validated local copy)")
    parser.add_argument("--loops", help="new loops output, from the `loops` subcommand")
    parser.add_argument("--scanlan-case", help="new Scanlan-Case output, from the `scanlan-case` subcommand")
    parser.add_argument("--analyze-output", help="new refined output, from the `analyze` subcommand")
    parser.add_argument("--analyze-raw", help="new raw output, from `analyze --raw`")
    args = parser.parse_args()

    strands_new = load(args.strands) if args.strands else load(HERE / "strands.pkl")
    danglings_new = load(args.danglings) if args.danglings else load(HERE / "danglings.pkl")

    if args.loops:
        compare_strands_danglings_loops(strands_new, danglings_new, load(args.loops))

    if args.scanlan_case:
        compare_scanlan_case(load(args.scanlan_case), None)
    elif args.analyze_output:
        refined = load(args.analyze_output)
        raw = load(args.analyze_raw) if args.analyze_raw else None
        if raw is not None:
            compare_strands_danglings_loops(raw.strands, raw.danglings, (raw.loops, raw.loop_orders))
        compare_scanlan_case(refined, raw)


if __name__ == "__main__":
    main()
