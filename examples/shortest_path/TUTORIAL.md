# Tutorial: shortest-path analysis (`--function shortest-path`)

Using this folder's denser network as the example. Also the tutorial for the
**default** full LAMMPS data file format; `examples/atrp/TUTORIAL.md` covers the
legacy pre-trimmed bond-list format and the full topology pipeline.

All commands assume you're in `src/` (paths are relative to there, matching
`src/README.md`); adjust if you run from elsewhere. Every command below was run
against the current code to confirm it still works exactly as shown.

## 0. What's in this folder

- `atrp_s3_sys3_0.data` - the topology: a full LAMMPS data file (`write_data` output),
  104,704 atoms, 114,029 bonds, 9,537 crosslinkers - denser than `examples/atrp`
  (2,932 crosslinkers), a good stress test for the crosslinker-pair search.
- `strands.pkl`, `danglings.pkl` - ground truth for `--function shortest-path`'s two
  searches, this tutorial's focus.
- `igloops.pkl`, `Scanlan_Case.pkl`, `CLUSTER_TESTING.md`, `compare_new_vs_old.py` -
  left over from testing the full topology pipeline on this denser network during
  development (loop detection here needed a bigger machine than a laptop - see
  `CLUSTER_TESTING.md` for that history). Not part of this tutorial; see
  `examples/atrp/TUTORIAL.md` if you want the full `--function topology` walkthrough.

## 1. What `--function shortest-path` computes, and why `--action` is required

This is **not** a generic "shortest path between any two nodes you name" tool - there
is no way to pass arbitrary source/target IDs. It only knows how to run exactly two
searches, and `--action` is how you pick which one - it isn't optional metadata, it's
the thing being computed, since it decides which nodes even count as a valid
source/target:

- `--action crosslinkers`: search between every pair of **crosslinker nodes**
  (degree >= 3). A path found this way is *called* a strand (one that doesn't pass
  through a third crosslinker along the way - if it did, that would actually be two
  strands, so those are excluded).
- `--action danglings`: search from every **crosslinker node** to every **chain-end
  node** (degree == 1). A path found this way is *called* a dangling chain - a defect,
  since one end isn't a load-bearing junction.

"Strand" and "dangling" aren't a separate concept layered on top of the search -
they're simply the names for what each action's result *is*, in polymer-network
terms. There's no third, generic output distinct from these two.

Both actions reuse the same multiprocessing search (`graph_analysis.crosslink`):
build the candidate pair set once, farm shortest-path search out across `--ncpus`
workers, cache the pair set under `cache-<sys-name>-*/` so a crash doesn't force
recomputing it. The output is always raw/unrefined - `--function shortest-path`
doesn't know about loops or the Scanlan-Case criterion at all; that's a separate
function (`--function topology`), which uses this same search internally as one of
its stages.

## 2. Run both searches

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
fast, since a path through a third crosslinker doesn't count); danglings is well under
a minute (far fewer chain ends than crosslinkers - 335 here, confirmed just now in
~43 seconds). Expect 17,934 strands and 335 danglings, matching this folder's shipped
files exactly.

### (Optional, unrelated) checking your output against the shipped files

`--function compare` has nothing to do with shortest-path specifically - it's a
generic utility that prints the entry count of any two pickle files, whatever
produced them. It's handy here purely as a quick sanity check, not part of the
shortest-path computation itself:

```bash
python3 main.py --function compare --file1 strands.pkl --file2 ../examples/shortest_path/strands.pkl
python3 main.py --function compare --file1 danglings.pkl --file2 ../examples/shortest_path/danglings.pkl
```

Confirmed `335 335` for danglings just now (two equal counts = good sign; it isn't a
byte-for-byte diff, just entry counts).

## 3. No `--format` flag needed here - that's the point

Unlike `examples/atrp` (which needs `--format bond-list`), this folder's topology file
is already the default full LAMMPS data file, so every command above omits
`--format` entirely. If you only have a bare, pre-trimmed bond table (just
`bond-id bond-type atom_i atom_j` rows), add `--format bond-list` - see
`examples/atrp/TUTORIAL.md`.

## 4. That's the whole function

`strands.pkl`/`danglings.pkl` are the complete output of `--function shortest-path` -
it doesn't do loop detection or Scanlan-Case refinement, by design (those are a
different function, `--function topology`, covered in `examples/atrp/TUTORIAL.md`).
If you want the fully refined result instead of just these two raw searches, run
`--function topology` directly on `atrp_s3_sys3_0.data` - it performs this same search
internally as one stage among several.
