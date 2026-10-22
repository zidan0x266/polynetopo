# Example 1: topology analysis (`--function topology`)

Walks through computing elastically effective strands, dangling chains, primary and
higher-order loops, and the Scanlan-Case crosslinker density refinement, for the ATRP
network in `examples/data/atrp_s1_sys1_0.data`. The other two examples use the same file.

All commands run from the repository root. Without `pip install -e .`, use
`python3 src/main.py` as below and put `src` on the path for your own scripts
(`PYTHONPATH=src`). On a cluster, load the Python module first (see the README's note on
import order).

## 0. What's here

- `examples/data/atrp_s1_sys1_0.data` - the network: a full LAMMPS data file
  (`write_data` output), 103,149 atoms and 105,869 bonds, after removing unreacted
  monomer at 96% conversion. Only its Bonds section is read by `--function topology`.
- `results.pkl` - the refined result (a `graph_analysis.RefinedStats`).
- `results_strands.pkl`, `results_danglings.pkl`, `results_loops.pkl` - the three raw
  stages, written by the same run so the refinement can be redone without the loop search.
- `run_log.txt` - the log and `/usr/bin/time -v` report of the cluster run that produced
  these files.

## 1. One command

```bash
python3 src/main.py --function topology \
    --input-file examples/data/atrp_s1_sys1_0.data \
    --output-file results.pkl
```

Strands, danglings, loops and the refinement come from one load of the file. The log
should end with:

```
Found 5156 strands, 395 danglings, 2721 loops
```

The strand and dangling searches take under a second; the loop search (igraph's minimum
cycle basis) is nearly all of the time and of the memory:

| | wall time | peak memory |
|---|---|---|
| cluster node, one core (`run_log.txt`) | 6 min 23 s | 10.4 GB |
| laptop | 2 min 24 s | 8.6 GB |

Add `--raw` to also write `results_raw.pkl`, a `graph_analysis.NetworkStats` holding the
same unrefined strands, danglings and loops in one object.

## 2. Check the numbers

```bash
PYTHONPATH=src python3 - <<'PY'
import igraph                      # first, see the README
import pickle
with open("results.pkl", "rb") as f:
    refined = pickle.load(f)
for name in ("effective_crosslinker_d4", "effective_crosslinker_d3", "all_crosslinker_d4",
             "all_crosslinker_d3", "correct_crosslinker_d4", "correct_crosslinker_d3",
             "switched_crosslinker_d4"):
    print(name, getattr(refined, name))
print("refined:", len(refined.strands), "strands,", len(refined.danglings), "danglings,", len(refined.loops), "loops")
PY
```

Expected:

```
effective_crosslinker_d4 2025
effective_crosslinker_d3 401
all_crosslinker_d4 2903
all_crosslinker_d3 29
correct_crosslinker_d4 878
correct_crosslinker_d3 -372
switched_crosslinker_d4 377
refined: 4650 strands, 398 danglings, 2721 loops
```

Before the refinement there are 5,156 strands and 395 danglings; of the 2,721 loops, 496
have exactly one crosslinker on them (primary loops).

## 3. One stage at a time

Each stage writes only its own file and never runs the refinement, because the files it
needs may not exist. The files `--stage full` wrote work as inputs here:

```bash
D=examples/data/atrp_s1_sys1_0.data
python3 src/main.py --function topology --stage strands   --input-file $D --output-file strands.pkl
python3 src/main.py --function topology --stage danglings --input-file $D --output-file danglings.pkl
python3 src/main.py --function topology --stage loops     --input-file $D --output-file igloops.pkl

# redo only the refinement, from the three files above
python3 src/main.py --function topology --stage scanlan-case --input-file $D \
    --strands-file strands.pkl --danglings-file danglings.pkl --loops-file igloops.pkl \
    --output-file Scanlan_Case.pkl
```

The refined result from these files is the same as `results.pkl`. This is the way to try
a different network quickly: run `strands` and `danglings` (seconds), and only run `loops`
(minutes, several GB) when you need it.

`python3 src/main.py --function compare --file1 strands.pkl --file2 results_strands.pkl`
prints the two entry counts as a quick check; comparing the contents needs a short script.

## 4. What each output means

- **strands / danglings** - flat lists of atom-ID paths. A strand's first and last IDs are
  crosslinkers (degree >= 3) and every atom between them has degree 2; the lower ID comes
  first. A dangling chain runs from a crosslinker to a free chain end (degree == 1).
  Several parallel chains between the same two crosslinkers are separate strands, and a
  direct bond between two crosslinkers is a strand of two atoms (360 here).
- **loops** - `[loops, crosslinker_counts]`: `loops[i]` is the member atoms of one
  independent smallest ring, `crosslinker_counts[i]` how many of them are crosslinkers.
  A count of 1 is a *primary loop*, a dead-end ring off a single junction - a defect, not a
  load path.
- **`results.pkl`** - the refined strands, danglings and loops after reclassifying
  crosslinkers that carry a primary loop or dangling chains as linear units, plus the
  resulting effective-vs-total crosslinker density (`effective_crosslinker_d4/d3`,
  `all_crosslinker_d4/d3`, ...). See `src/graph_analysis/scanlan_case.py`.

## 5. Notes

- **Node IDs are the raw atom IDs of the data file** in every output, so a loop's atoms
  can be matched to strand and dangling endpoints directly.
- **Loops depend on edge order.** A minimum cycle basis is not always unique, and igraph's
  choice among equally short cycles depends on the order the bonds were added, which is the
  order of the Bonds section. The total loop count and the primary-loop count do not change,
  but a few individual loops can.
- **Older outputs.** Earlier releases found 5,138 strands here, not 5,156: their search
  missed strands when three or more chains join the same two crosslinkers (see the main
  README). The crosslinker densities above are the same either way.
- **Legacy inputs.** `--format bond-list` and `--format rmc` read a bare bond table or an
  RMC-style edge list; there is no example file for them, since the LAMMPS reader makes
  them unnecessary here.

See the main [README](../../README.md) and the topology papers cited there for the
scientific background (the Scanlan-Case criterion, what makes a crosslinker "effective").
