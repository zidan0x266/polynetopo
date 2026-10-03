# Tutorial: spatial betweenness analysis (`--function betweenness`)

Walks through computing GEBC (geodesic edge betweenness centrality), OBC
(orientation-weighted) and DBC (distance-weighted) for every bond in one pass, using
this folder's two network states as the example.

## 0. What's in this folder

- `atrp_s1_sys1_0.data` - the **undeformed** state: 103,149 atoms, 105,869 bonds, a
  cubic box (all three box lengths equal, ~49.397).
- `atrp_s1_sys1_42.data` - a **deformed** state of the same network: same atom/bond
  counts, a non-cubic box (stretched along one axis). Good for seeing how OBC/DBC
  respond to a loading direction on a real deformed geometry, not just the
  undeformed reference.

Neither file has precomputed centrality output shipped in this folder (unlike
`examples/atrp`/`examples/shortest_path`) - this tutorial generates it.

## 1. Build the in-house igraph (once per machine)

`--function betweenness` needs `Graph.edge_betweenness_spatial`, a method the stock
PyPI `igraph` package doesn't have. Build it once:

```bash
git submodule update --init vendor/igraph vendor/python-igraph
echo "1.0.0" > vendor/igraph/IGRAPH_VERSION   # the fork has no git tags; see vendor/README.md
cd vendor && bash python-igraph/scripts/build_spatial.sh && cd ..
vendor/.spatial-venv/bin/python -m pip install --no-build-isolation --no-deps networkx tqdm
vendor/.spatial-venv/bin/python -m pip install -e .
```

Takes ~15-20 minutes; full details (including the exact CMake error you'd hit without
the `IGRAPH_VERSION` line) are in `../../vendor/README.md`. Confirm it worked:

```bash
vendor/.spatial-venv/bin/python -c 'import igraph as ig; print(hasattr(ig.Graph, "edge_betweenness_spatial"))'
# must print True
```

Every command below uses `vendor/.spatial-venv/bin/python`, never your regular
environment's `python3` - stock `igraph` doesn't have the native method, and calling
`--function betweenness` with it raises a clear `RuntimeError` rather than silently
doing the wrong thing.

## 2. Inspect before computing anything

```bash
vendor/.spatial-venv/bin/python src/main.py --function betweenness \
    --input-file examples/betweenness/atrp_s1_sys1_0.data --inspect
```

Parses and validates the file, prints its topology and *actual* box lengths as JSON,
and stops - no centrality calculated. Confirms the undeformed file's box really is
cubic (49.397 x 49.397 x 49.397) before you commit to a multi-minute run. Try the
deformed file too - its box is not cubic, which matters for OBC/DBC.

## 3. Fast sanity check on a small sample

```bash
vendor/.spatial-venv/bin/python src/main.py --function betweenness \
    --input-file examples/betweenness/atrp_s1_sys1_0.data \
    --sample-vertices 256 --validate-stock --output-file /tmp/sample.pkl --overwrite
```

Seconds, not minutes - but a genuinely different, smaller induced graph, not an
estimate of the full network. `--validate-stock` runs a separate stock
`igraph.Graph.edge_betweenness()` pass and fails loudly if native GEBC disagrees
(`rtol=1e-12, atol=1e-8`); on this network it agrees with
`max_absolute_difference: 0.0`.

## 4. Full calculation, with first-run validation

```bash
vendor/.spatial-venv/bin/python src/main.py --function betweenness \
    --input-file examples/betweenness/atrp_s1_sys1_0.data --direction x --validate-stock
```

Validating roughly doubles runtime (a second full traversal) - worth it once per new
network/machine combination; on this 103,149-atom/105,869-bond network, expect on the
order of ten minutes total with `--validate-stock`, proportionally less without it.
Writes `atrp_s1_sys1_0.x.spatial.pkl` next to the input (omit `--output-file` to get
this default, or set one explicitly).

For the deformed file, pass its actual loading direction and skip re-validating
(already proven on this machine in step 3/4 above):

```bash
vendor/.spatial-venv/bin/python src/main.py --function betweenness \
    --input-file examples/betweenness/atrp_s1_sys1_42.data --direction x
```

## 5. Read the result

Output is a pickled `graph_analysis.betweenness.CentralityResult` - a plain pickle of
a dataclass, like every other output in this package, not a NumPy archive:

```python
import pickle
with open("examples/betweenness/atrp_s1_sys1_0.x.spatial.pkl", "rb") as f:
    result = pickle.load(f)  # needs graph_analysis importable: pip install -e .

bond_ids = result.network.bond_ids        # aligned with every array below
raw_gebc, raw_obc, raw_dbc = result.raw_gebc, result.raw_obc, result.raw_dbc
obc, dbc = result.obc, result.dbc         # pair-normalized properties, not stored separately
print(result.summary())                   # timing, validation, box, neff, ...
```

`raw_dbc` is **always** the distance-only sum - there's no automatic `Lx` or
bond-length division. Apply any reference length yourself, as many times as you like,
without recomputing:

```python
dbc_undeformed_ref = result.dbc_for_length(49.397)       # raw_dbc / 49.397, unnormalized
dbc_undeformed_ref_norm = result.dbc_for_length(49.397, normalized=True)  # + pair normalization
```

`--dbc-length L` at the CLI does exactly this, applied once and stored in
`result.dbc`/the saved file, so you don't have to recompute it yourself if one
reference length is all you need. See `vendor/python-igraph/SPATIAL_BETWEENNESS.md`
and `vendor/python-igraph/analysis/README.md` for the full scientific-convention
reference and the extensive validation already recorded there (that documentation's
examples predate this pickle-based output and show `.npz` field names - they map
directly onto this pickle's attributes, e.g. `raw_GEBC_igraph`/`OBC`/`DBC` there are
`result.raw_gebc`/`obc`/`dbc` here).

## 6. Compare the two states

Since both files are the same network at two points along a deformation, a quick way
to see the effect of stretching is to compare their GEBC/OBC/DBC distributions
directly (same `bond_ids` in both, since bonds aren't created/destroyed between these
two states):

```python
import pickle
import numpy as np

with open("examples/betweenness/atrp_s1_sys1_0.x.spatial.pkl", "rb") as f:
    undeformed = pickle.load(f)
with open("examples/betweenness/atrp_s1_sys1_42.x.spatial.pkl", "rb") as f:
    deformed = pickle.load(f)

assert np.array_equal(undeformed.network.bond_ids, deformed.network.bond_ids)
print("mean OBC, undeformed vs deformed:", undeformed.obc.mean(), deformed.obc.mean())
print("mean DBC, undeformed vs deformed:", undeformed.dbc.mean(), deformed.dbc.mean())
```

This is exactly the kind of comparison the fracture-prediction paper in
`../../references/` is built on: OBC from the *initial* state predicts final fracture
locations best, while DBC correlates more strongly with fractures later in stretching
- see that paper for the full picture.
