# Example 2: spatial betweenness (`--function betweenness`)

Computes, for every bond of the network in `examples/data/atrp_s1_sys1_0.data`, the
geodesic edge betweenness (GEBC) and its two spatially weighted variants, the
orientation-weighted OBC and the distance-weighted DBC (Bai et al., Macromolecules 2026),
for the loading axis x. All three come out of one pass.

All commands run from the repository root. This function needs the in-house igraph build
(`Graph.edge_betweenness_spatial`); check it with
`python3 -c "import igraph; print(hasattr(igraph.Graph, 'edge_betweenness_spatial'))"`,
which must print `True`. On a cluster, load the Python module first (see the README's note
on import order).

## 0. What's here

- `atrp_s1_sys1_0.x.spatial.pkl` - the result, a pickled `CentralityResult` (9 MB).
- `summary.json` - the run's metadata and the comparison with stock igraph.
- `run_log.txt` - the log and `/usr/bin/time -v` report of the cluster run that produced
  them.

## 1. Inspect the input first

```bash
python3 src/main.py --function betweenness \
    --input-file examples/data/atrp_s1_sys1_0.data --inspect
```

This parses and validates the file, prints its topology and box as JSON, and computes
nothing. Expect `"nodes": 103149`, `"edges": 105869`, `"components": 1`, and box lengths of
49.3969 on each axis. Always look at the box lengths before a run on a new network, since
OBC and DBC depend on them. The box is taken as periodic along all three axes.

## 2. A quick trial

```bash
python3 src/main.py --function betweenness \
    --input-file examples/data/atrp_s1_sys1_0.data --direction x \
    --sample-vertices 256 --validate-stock --output-file sample.pkl
```

This runs on a 256-atom neighbourhood (257 bonds), a different, smaller graph and not an
estimate of the full result, and takes seconds. Expect `"allclose": true` and, in the
result, a largest raw GEBC of 15159.0, OBC of 7173.69 and DBC of 111777.42 (values from a
laptop run).

## 3. The full run

```bash
python3 src/main.py --function betweenness \
    --input-file examples/data/atrp_s1_sys1_0.data --direction x --validate-stock \
    --output-file atrp_s1_sys1_0.x.spatial.pkl
```

The calculation is single-threaded, so more cores don't speed it up, and it needs little
memory:

| | spatial pass | `--validate-stock` pass | total | peak memory |
|---|---|---|---|---|
| cluster node (`run_log.txt`) | 43 min | 22 min | 1 h 06 min | 119 MiB |
| laptop | 24 min | 9 min | 34 min | 94 MiB |

The cluster's spatial pass shared the node with two other jobs for its first 17 minutes,
which made that part about 2.5 times slower. Its steady pace afterwards, 89 to 107 s per
5%, puts a pass on a quiet node at about 30 minutes.

`--validate-stock` also computes the GEBC with stock igraph's `edge_betweenness()` and
stops with an error if the two disagree. It is a second full traversal, worth doing once
for each new network or machine and then leaving out. Here it agreed exactly:
`"max_absolute_difference": 0.0` and `"allclose": true` in `summary.json`.

## 4. Read the result

```bash
PYTHONPATH=src python3 - <<'PY'
import igraph                      # first, see the README
import pickle
with open("examples/example_2/atrp_s1_sys1_0.x.spatial.pkl", "rb") as f:
    result = pickle.load(f)
print(result.raw_gebc.shape, result.neff, result.raw_gebc.max(), result.raw_obc.max(), result.raw_dbc.max())
PY
```

Expected: `(105869,) 103149`, then a largest raw GEBC of 1.07898e+08, OBC of 6.13041e+07 and
DBC of 2.77176e+09. The largest GEBC is on bond 69986.

- `result.raw_gebc`, `raw_obc` and `raw_dbc` have one value per bond, in the order of
  `result.network.bond_ids`.
- `result.obc` and `result.dbc` are the same values divided by the number of atom pairs
  (`normalization_factor = 2 / (neff (neff - 1))`, here 1.88e-10). They are computed when you
  ask for them, not stored.
- `result.dbc_for_length(L)` divides the raw DBC by a reference length `L` of your choice,
  for example the box length, and `dbc_for_length(L, normalized=True)` also normalizes it.
  `--dbc-length L` at the command line changes only the `dbc` property, never `raw_dbc`.

## 5. What the three numbers are

For an edge e, a sum over atom pairs (s, t) of the fraction of shortest s-t paths that
pass through e, times a weight w_st:

| | w_st |
|---|---|
| GEBC | 1 |
| OBC | cos of the angle between the s-t displacement and the loading axis, \|dx\| / \|r\| |
| DBC | the distance \|r\| between s and t (the paper divides by the box length L_x; the raw value here does not) |

Each unordered pair is counted once, shortest paths are unweighted (every bond counts 1),
and the s-t displacement is the minimum image under the periodic box.

## 6. How it was checked

- **GEBC against stock igraph**, on the full network, exactly as above.
- **All three against a brute-force implementation** of these definitions, written without
  using the package, on random periodic networks of 90 to 100 atoms along each axis:
  relative differences near 1e-15.
- **The cluster result against an independent laptop run** of the same network: the GEBC is
  bit-for-bit identical and the OBC and DBC agree to within 6e-16 relative.

## 7. Notes

- A pickle written by an earlier release of the package, from when `betweenness` was a
  folder, no longer loads.
- To run the same analysis for the other axes, change `--direction`; each axis is a separate
  full run.
