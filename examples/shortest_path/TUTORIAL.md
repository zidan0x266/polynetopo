# Tutorial: shortest-path analysis (`--function shortest-path`)

Using this folder's denser network as the example. Also the tutorial for the
**default** full LAMMPS data file format; `examples/atrp/TUTORIAL.md` covers the
legacy pre-trimmed bond-list format and the full topology pipeline in more depth.

All commands assume you're in `src/` (paths are relative to there, matching
`src/README.md`); adjust if you run from elsewhere. Every command below was run
against the current code to confirm it still works exactly as shown.

## 0. What's in this folder

- `atrp_s3_sys3_0.data` - the topology: a full LAMMPS data file (`write_data` output),
  104,704 atoms, 114,029 bonds, 9,537 crosslinkers - denser than `examples/atrp`
  (2,932 crosslinkers), a good stress test for the crosslinker-pair search.
- `strands.pkl`, `danglings.pkl`, `igloops.pkl`, `Scanlan_Case.pkl` - ground truth for
  every stage, not just shortest-path, since this folder doubles as the full-pipeline
  validation example (see `CLUSTER_TESTING.md` for how the loop-detection stage was
  validated on a bigger machine than this one).
- `compare_new_vs_old.py` - a ready-made script that diffs a fresh run's
  strands/danglings/loops/Scanlan-Case against the files above.

## 1. If you want the complete, SC-refined result: one command

If what you actually want is "the shortest-path calculation" in the broader sense -
a finished, Scanlan-Case-refined set of strands/danglings/loops - that's
`--function topology`, not `--function shortest-path`. It runs the crosslinker and
dangling-chain searches below internally, as one of its stages, then refines the
result automatically - no separate strand/dangling/loop files to juggle:

```bash
python3 main.py --function topology \
    --input-file ../examples/shortest_path/atrp_s3_sys3_0.data \
    --output-file results.pkl --raw --ncpus 4 --max-chunk-size 10000
```

`results.pkl` is a `graph_analysis.RefinedStats`: refined `strands`, `danglings`,
`loops`, `loop_orders`, all after SC refinement, plus the crosslinker-density
numbers - always, regardless of `--raw`. `--raw` additionally writes
`results_raw.pkl` (the pre-refinement snapshot, for reference/debugging only, not
required to get the three refined components above). See
`examples/atrp/TUTORIAL.md` for the full walkthrough of this pipeline, including
what each field means and the `--stage loops`/`scanlan-case` resumable variants.

The rest of this tutorial covers `--function shortest-path` specifically: the
lower-level crosslinker-to-crosslinker / crosslinker-to-chain-end search that
`--function topology` runs internally, exposed as its own standalone tool for when
you want just the raw strands or danglings, or want to run that one stage alone
(e.g. to checkpoint it separately on a memory-constrained machine).

## 2. What `--function shortest-path` computes on its own

Two independent searches, selected by `--action`:

- `--action crosslinkers`: every shortest simple path between two crosslinkers
  (degree >= 3) - an *elastically effective strand* if the path doesn't pass through
  a third crosslinker along the way.
- `--action danglings`: every shortest simple path from a crosslinker to a free chain
  end (degree == 1) - a *dangling chain*, a defect that doesn't bear load.

Both reuse the same multiprocessing search (`graph_analysis.crosslink`): build the
candidate pair set once, farm shortest-path search out across `--ncpus` workers, cache
the pair set under `cache-<sys-name>-*/` so a crash doesn't force recomputing it. This
output is always raw/unrefined - Scanlan-Case hasn't run yet at this stage.

## 3. Run both searches

```bash
# Strands: crosslinker-to-crosslinker
python3 main.py --function shortest-path --sys-name strands \
    --input-file ../examples/shortest_path/atrp_s3_sys3_0.data \
    --output-file strands.pkl --action crosslinkers --ncpus 4 --max-chunk-size 10000
rm -r cache-*

# Dangling chains: crosslinker-to-chain-end
python3 main.py --function shortest-path --sys-name danglings \
    --input-file ../examples/shortest_path/atrp_s3_sys3_0.data \
    --output-file danglings.pkl --action danglings --ncpus 4 --max-chunk-size 10000
rm -r cache-*
```

On this network: ~1 minute for strands (4.3M candidate crosslinker pairs - most fail
fast, since a path through a third crosslinker doesn't count), well under a minute for
danglings (far fewer chain ends than crosslinkers - 335 here, confirmed just now in
~43 seconds). Expect 17,934 strands and 335 danglings, matching this folder's shipped
files exactly:

```bash
python3 main.py --function compare --file1 strands.pkl --file2 ../examples/shortest_path/strands.pkl
python3 main.py --function compare --file1 danglings.pkl --file2 ../examples/shortest_path/danglings.pkl
```

(`compare` only prints the two counts - confirmed `335 335` just now for danglings;
for a true byte-for-byte diff, see `compare_new_vs_old.py` in this folder, which also
handles path-direction normalization and the known Scanlan-Case bug-fix delta.)

## 4. No `--format` flag needed here - that's the point

Unlike `examples/atrp` (which needs `--format bond-list`), this folder's topology file
is already the default full LAMMPS data file, so every command above omits
`--format` entirely. If you only have a bare, pre-trimmed bond table (just
`bond-id bond-type atom_i atom_j` rows), add `--format bond-list` - see
`examples/atrp/TUTORIAL.md`.

## 5. Finishing the pipeline from these two files

If you ran step 3 (standalone) rather than step 1 (one-shot), you still need loops +
Scanlan-Case to get a refined result:

```bash
python3 main.py --function topology --stage loops \
    --input-file ../examples/shortest_path/atrp_s3_sys3_0.data --output-file igloops.pkl

python3 main.py --function topology --stage scanlan-case \
    --input-file ../examples/shortest_path/atrp_s3_sys3_0.data \
    --strands-file strands.pkl --danglings-file danglings.pkl --loops-file igloops.pkl \
    --output-file Scanlan_Case.pkl --raw
```

This network's loop-detection stage is expensive enough that it needed a separate,
larger-memory machine to complete - see `CLUSTER_TESTING.md` for that story and the
exact commands used, and the main [README](../../README.md)'s reproducibility notes
for what that run revealed about loop detection in general (vertex-count padding is
safe; edge order is what actually affects which minimum cycle basis you get).
