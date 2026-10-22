# Tutorial: topology analysis (`--function topology`)

Walks through computing elastically effective strands, dangling chains, primary/
higher-order loops, and Scanlan-Case crosslinker density refinement, using this
folder's ATRP network as the example. This is also the tutorial for the **legacy**
pre-trimmed bond-list format - `examples/shortest_path/TUTORIAL.md` uses the default
full LAMMPS data file format instead.

All commands assume you're in `src/` (paths are relative to there, matching
`src/README.md`); adjust if you run from elsewhere.

## 0. What's in this folder

- `edges_0.bds` - the topology: a bare bond table (`bond-id bond-type atom_i atom_j`
  per line), 103,149 atoms after removing unreacted monomer at 96% conversion. This
  is the **legacy** format (`--format bond-list`); it was manually trimmed from a full
  LAMMPS data file before this package could read one directly.
- `strands.pkl`, `danglings.pkl`, `igloops.pkl`, `Scanlan_Case.pkl` - the expected
  output of every command below, already computed - use them to check your own run
  matches, no need to treat your first run as a leap of faith.

## 1. One-shot: the whole pipeline at once

```bash
python3 main.py --function topology \
    --input-file ../examples/atrp/edges_0.bds --format bond-list \
    --output-file results.pkl --raw
```

Strands, danglings, loops and their Scanlan-Case refinement, from one topology load,
no intermediate files. Strands and danglings take a few seconds on this network
(one pass over the bonds, see `src/graph_analysis/crosslink.py`); loop detection
(the minimum cycle basis) is the slow stage.

Check it against the shipped results:

```bash
python3 -c "
import pickle
with open('results.pkl', 'rb') as f: new = pickle.load(f)
with open('../examples/atrp/Scanlan_Case.pkl', 'rb') as f: shipped = pickle.load(f)
fields = ['effective_crosslinker_d4', 'effective_crosslinker_d3', 'all_crosslinker_d4',
          'all_crosslinker_d3', 'correct_crosslinker_d4', 'correct_crosslinker_d3',
          'switched_crosslinker_d4']
print('crosslinker-density stats match:', all(getattr(new, k) == getattr(shipped, k) for k in fields))
"
```

Expect `True` - `[2025, 401, 2903, 29, 878, -372, 377]`. The shipped
`strands.pkl` has 5,156 strands (earlier releases found 5,138 and missed some
parallel chains - see the main README's note on earlier strand counts); the density
stats above are unchanged by that. Danglings (395) and loops are unchanged.

## 2. Step-by-step: the same result, one stage at a time

Useful for inspecting an intermediate stage, or resuming after a crash - each stage
writes its own file.

```bash
# Strands: chains between two crosslinkers (degree >= 3)
python3 main.py --function shortest-path \
    --input-file ../examples/atrp/edges_0.bds --format bond-list \
    --output-file strands.pkl --action crosslinkers

# Dangling chains: chains from a crosslinker to a chain end (degree == 1)
python3 main.py --function shortest-path \
    --input-file ../examples/atrp/edges_0.bds --format bond-list \
    --output-file danglings.pkl --action danglings

# Minimum cycle basis (loops), with per-loop crosslinker counts
python3 main.py --function topology --stage loops \
    --input-file ../examples/atrp/edges_0.bds --format bond-list \
    --output-file igloops.pkl --nvertices 106250

# Scanlan-Case refinement of the three files above
python3 main.py --function topology --stage scanlan-case \
    --input-file ../examples/atrp/edges_0.bds --format bond-list \
    --strands-file strands.pkl --danglings-file danglings.pkl --loops-file igloops.pkl \
    --output-file Scanlan_Case.pkl --raw
```

Each output file should match this folder's shipped copy exactly
(`python3 main.py --function compare --file1 strands.pkl --file2 ../examples/atrp/strands.pkl`
prints two equal counts as a quick check; byte-for-byte comparison needs a small
script like the one in step 1).

`--nvertices 106250` above matches the number this network had *before* removing
unreacted monomer (106,250 beads, down to the 103,149 in `edges_0.bds` after
deletion) - a holdover from an earlier pipeline stage. It's safe either way: omitting
`--nvertices` (letting it infer from the topology) gives an identical result, since
vertex padding doesn't affect loop detection, only edge order does (both the full
data file and this bond-list format preserve the file's original line order).

## 3. What each output means

- **`strands.pkl`** / **`danglings.pkl`**: flat lists of atom-ID paths. A strand's
  first and last IDs are crosslinkers (degree >= 3); a dangling chain's last ID is a
  free chain end (degree == 1).
- **`igloops.pkl`**: `[loops, crosslinker_counts]` - `loops[i]` is the member-atom set
  of one independent smallest ring; `crosslinker_counts[i]` is how many of those atoms
  are crosslinkers. A count of 1 is a *primary loop* (a dead-end ring off a single
  junction - a defect, not a load path).
- **`Scanlan_Case.pkl`** (now a `graph_analysis.RefinedStats`, not the old nested-list
  format): refined strands/danglings/loops after reclassifying crosslinkers that
  carry a primary loop or dangling chain(s) as linear units, plus the resulting
  effective-vs-total crosslinker density (`effective_crosslinker_d4/d3`,
  `all_crosslinker_d4/d3`, ...). `--raw` additionally saves `Scanlan_Case_raw.pkl` (a
  `NetworkStats`): the pre-refinement strands/danglings/loops, for sanity-checking a
  configuration without redoing the analysis.

See the main [README](../../README.md)'s reproducibility notes and the two topology
papers in `../../references/` for the full scientific background (the Scanlan-Case
criterion, what makes a crosslinker "effective").
