# Example 3: shortest paths across the periodic boundary (`--function shortest-path`)

The method of Yu & Jackson, "Shortest Paths Govern Bond Rupture in Thermoset Networks"
(Macromolecules 2025), run on `examples/data/atrp_s1_sys1_0.data` along x. The README's
section on this function has the definitions; in short, for each period p = 1..M, every
atom's shortest path to its own image p box lengths away is found on a graph of p + 1
periodic images, and its length divided by p is the normalized length D^g (bonds per box
length). The shortest one predicts where the first bond breaks and, through D^g_min / L0,
the strain at which it breaks.

All commands run from the repository root. On a cluster, load the Python module first
(see the README's note on import order).

## 0. What's here

- `atrp_s1_sys1_0.x.sp.pkl` - the result, a pickled `ShortestPathResult` (20 MB).
- `summary.json` - the numbers the command prints.
- `shortest_path_lengths.csv` - the 3,321 normalized lengths D^g, one per line, so the large
  pickle is optional.
- `run_log.txt` - the log and `/usr/bin/time -v` report of the cluster run that produced
  them.

## 1. A quick trial

```bash
python3 src/main.py --function shortest-path \
    --input-file examples/data/atrp_s1_sys1_0.data --direction x \
    --slice-max 1.0 --max-period 2 --ncpus 4 --output-file trial.sp.pkl
```

`--slice-max 1.0` starts paths only from atoms within 1.0 of the box's low face (2,069 of
103,149 atoms), which makes it a small test and not a result for the network. Expect:

```
period 1: 706 valid paths from 2069 sources
period 2: 456 valid paths from 2069 sources
Found 298 shortest paths
```

and in the printed summary `"min_length": 155.0`, `"max_length": 419.0`.

## 2. The full run

```bash
python3 src/main.py --function shortest-path \
    --input-file examples/data/atrp_s1_sys1_0.data --direction x --ncpus 8 \
    --output-file atrp_s1_sys1_0.x.sp.pkl
```

With the default `--max-period 8` the search stops at period 6, the first with no valid
path. The log shows how many of the 103,149 atoms found a valid path at each period:

| period p | valid paths | time (34 workers) |
|---|---|---|
| 1 | 34,534 | about 1 min |
| 2 | 22,353 | 2 min 05 s |
| 3 | 10,535 | 2 min 38 s |
| 4 | 1,997 | 3 min 09 s |
| 5 | 66 | 3 min 43 s |
| 6 | 0 | 4 min 16 s |

Removing redundant paths then leaves **3,321 paths**, and the summary is `"min_length": 155.0`,
`"max_length": 419.0`, `"mean_length": 218.80`, `"median_length": 211.0`.

| | wall time | peak memory of the largest process |
|---|---|---|
| cluster node, 34 workers (`run_log.txt`) | 17 min | 0.79 GB |
| laptop, 8 workers (earlier version of the code) | 106 min | 5.8 GB |

Each worker holds one copy of the multi-image graph, so memory grows with `--ncpus` and with
the number of periods. The time per period grows with p because the paths and graphs do.

## 3. Read the result

```bash
PYTHONPATH=src python3 - <<'PY'
import igraph                      # first, see the README
import pickle
with open("examples/example_3/atrp_s1_sys1_0.x.sp.pkl", "rb") as f:
    sp = pickle.load(f)
print(len(sp.paths), sp.min_length, sp.min_length_by_period.tolist())
print(sp.paths[0][:3].tolist(), "...", sp.paths[0][-1].tolist())
PY
```

Expected: `3321 155.0`, then the smallest D^g per period
`[[1.0, 155.0], [2.0, 173.0], [3.0, 168.33], [4.0, 172.0], [5.0, 174.0]]`, then
`[[6, 0], [28, 0], [219, 0]] ... [6, 1]` for the first path.

- `sp.lengths[k]` is the D^g of `sp.paths[k]`.
- A path is an array of `(atom ID, image)` rows from the atom in image 0 to the same atom in
  image p. This path starts and ends on atom 6, at image 0 and image 1, and has 346 bonds.
- Images are counted from image 0, the box as given. An image of -1, or above p, is one of
  the buffer images a path passed through; here the paths reach images -1 to 5.
- `sp.min_length_by_period` is computed before redundant paths are removed.

## 4. The options

- `--max-period M` (default 8): periods 1..M. Here nothing is found past period 5, so 8 and 6
  give the same result.
- `--buffer-images N` (default 1): extra images on each side of the p + 1, so a path can leave
  the stack between its ends. `0` keeps it inside, and the paths can differ.
- `--min-sp-only`: also apply the paper's optional Case 3 and keep only the shortest path of
  each start atom. The 3,321 paths here start from 2,526 different atoms, so this would keep
  2,526 of them, with the same minimum of 155.0 (counted by applying the same step to this
  result, not from a separate run). It is off by default, so an atom can have a path at
  several periods.
- `--keep-redundant`: skip the removal of redundant paths (the paper's Cases 1 and 2).
- `--slice-max L`, `--ncpus N`, `--checkpoint-file F`: start paths from a thin slice, use
  worker processes, and save progress after each period so an interrupted run resumes.

## 5. How it was checked

On this network the lengths and the paths, row by row, are identical to those of the original
code, run with the same buffer image, on the cluster: the same 3,321 lengths, from 155 to 419.
On random periodic networks, the results equal those of the original code for
buffer images 0 to 2, along all three axes, with one and with several workers, and for Case 3.

## 6. Notes

- The input may not contain an atom pair bonded twice, or an atom bonded to itself: a path
  lists atoms, so it could not say which of two parallel bonds it used. This raises an error
  naming the pair. The reaction code that produced this network does not create such bonds.
- Atom IDs are the raw IDs of the data file.
- Run `--direction y` and `--direction z` for the other axes; each is a separate run.
