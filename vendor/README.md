# In-house igraph (spatial betweenness)

`--function betweenness` (GEBC/OBC/DBC) needs the in-house igraph build here,
not the stock PyPI `igraph` package used by `--function topology`/`shortest-path`.
These are git submodules (two sibling repos, forked at the official `1.0.0` tag
on a `spatial-betweenness-1.0.0` branch that adds one native method,
`Graph.edge_betweenness_spatial`), not vendored copies - cloning this repo does
not pull their content until you run `git submodule update --init`.

## Build (once per machine - binaries aren't portable across machines/OSes)

```bash
git submodule update --init vendor/igraph vendor/python-igraph

# This fork has no git tags (lost when it was forked), and igraph's CMake build
# needs one to determine its version string. Bypass with the known version:
echo "1.0.0" > vendor/igraph/IGRAPH_VERSION

cd vendor
bash python-igraph/scripts/build_spatial.sh
```

This creates `.spatial-venv` (a dedicated virtualenv) and `.spatial-build`
(C build artifacts), both gitignored. It builds the C core, runs its test
suite, builds the Python extension against it, and runs the Python test suite -
expect ~15-20 minutes. Then add this project's own dependencies to the same
venv so it becomes the one environment for every `--function`:

```bash
.spatial-venv/bin/python -m pip install --no-build-isolation --no-deps networkx tqdm
.spatial-venv/bin/python -m pip install -e ..   # the graph_analysis package itself
```

Use `vendor/.spatial-venv/bin/python src/main.py --function betweenness ...` (or
any other `--function`) from the repo root. See
`vendor/python-igraph/SPATIAL_BETWEENNESS.md` for the full build reference,
scientific conventions (raw DBC is distance-only; apply a reference length
yourself - see `CentralityResult.dbc_for_length`), and the extensive validation
already recorded there against stock igraph and independent path enumeration.
