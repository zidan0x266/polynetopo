# Polymer network topology, shortest-path and betweenness analysis

A single CLI (`graph_analysis.cli`, exposed here as `main.py`) covering the
crosslinked-network analysis pipeline: elastically effective strands, dangling
chains, primary/higher-order loops, Scanlan-Case crosslinker density refinement
(used for the ATRP/FRP network papers), and native GEBC/OBC/DBC spatial edge
betweenness (used for the fracture-prediction paper) - see `../references`.

One flag, `--function {topology, shortest-path, betweenness, compare}` (default
`topology`), selects what to run; run `python3 main.py --help` for the full option
list, grouped by which function each belongs to.

## Topology formats

`--format` (default `lammps`) selects how `--input-file` is read, for
`--function topology`/`shortest-path`:

- `lammps` (default) - a full LAMMPS data file (`write_data` output, with its
  Masses/Atoms/Velocities/Bonds/... sections); only the Bonds section is read.
- `bond-list` - **legacy**: a bare bond table already trimmed down to just
  `bond-id bond-type atom_i atom_j` per line, nothing else. Kept for topology files
  saved this way before this package could read a full data file directly, and as a
  starting point to adapt for a bond list from a different tool (e.g. a reactive
  Monte Carlo / 3DBFM code) if its column layout differs.
- `rmc` - **legacy**: an RMC-style whitespace edge list.

`--function betweenness` always reads a full LAMMPS data file directly (see below).

## `--function topology` (one-shot, recommended for a new project)

`--stage full` (the default) runs every stage - strands, danglings, loops,
Scanlan-Case refinement - from a single topology load, with no intermediate files:

```bash
python3 main.py --function topology \
    --input-file ../examples/shortest_path/atrp_s3_sys3_0.data \
    --output-file results.pkl --raw --ncpus 4 --max-chunk-size 10000
```

Writes `results.pkl` (a `graph_analysis.RefinedStats`: the Scanlan-Case-refined
strands/danglings/loops plus the crosslinker density it produced). Add `--raw` to
also write `results_raw.pkl` (a `graph_analysis.NetworkStats`: the pre-refinement
strands/danglings/loops, for sanity-checking a specific configuration without
redoing the analysis). `--nvertices` is optional - if omitted it's inferred from the
topology (see the reproducibility note in `graph_analysis/loops.py` before relying on
this to reproduce an old run bit-for-bit: vertex count doesn't affect the result, but
edge *order* does, which is why every topology loader preserves the file's original
line order rather than going through `networkx`'s own edge iteration).

### `--stage loops`/`scanlan-case` (step-by-step, resumable)

Useful for inspecting an intermediate stage, running a memory-heavy stage on a bigger
machine separately, or resuming after a crash - each stage caches its crosslinker
search matrix under `cache-<sys-name>-*/`. The `examples/atrp` dataset only has a
pre-trimmed bond list, so it needs `--format bond-list`:

```bash
# Strands: shortest paths between crosslinkers (degree >= 3)
python3 main.py --function shortest-path --sys-name strands --input-file ../examples/atrp/edges_0.bds \
    --output-file strands.pkl --action crosslinkers --max-chunk-size 10000 --ncpus 4 --format bond-list
rm -r cache-*

# Dangling chains: shortest paths from a crosslinker to a chain end (degree == 1)
python3 main.py --function shortest-path --sys-name danglings --input-file ../examples/atrp/edges_0.bds \
    --output-file danglings.pkl --action danglings --ncpus 4 --max-chunk-size 10000 --format bond-list
rm -r cache-*

# Minimum cycle basis (loops), with per-loop crosslinker counts
python3 main.py --function topology --stage loops --input-file ../examples/atrp/edges_0.bds \
    --output-file igloops.pkl --format bond-list --nvertices 106250

# Scanlan-Case refinement of the strands/danglings/loops above (also takes --raw)
python3 main.py --function topology --stage scanlan-case --input-file ../examples/atrp/edges_0.bds \
    --strands-file strands.pkl --danglings-file danglings.pkl --loops-file igloops.pkl \
    --output-file Scanlan_Case.pkl --format bond-list

# Compare the entry counts of two pickled result files
python3 main.py --function compare --file1 strands.pkl --file2 danglings.pkl
```

`examples/shortest_path` has the full LAMMPS data file, so the same commands work
with `--input-file ../examples/shortest_path/atrp_s3_sys3_0.data` and no `--format`
flag needed (it's the default).

`Scanlan_Case.pkl`/`results.pkl` now pickle a `graph_analysis.RefinedStats` dataclass
rather than the old positional nested-list structure; call `.to_legacy()` on the
`ScanlanCaseAnalyzer.run()` return value (before pickling) to recover the previous
`[[refined...], [raw...], [counts...]]` shape if needed for older downstream scripts.

## `--function betweenness`

Native GEBC/OBC/DBC edge betweenness, computed together in one pass. Needs the
in-house igraph build in `../vendor/` (see `../vendor/README.md`) - run with that
build's own Python, not your regular environment's:

```bash
../vendor/.spatial-venv/bin/python main.py --function betweenness \
    --input-file ../examples/betweenness/atrp_s1_sys1_0.data --direction x --validate-stock
```

- `--input-file` is always a full LAMMPS data file (`Atoms # full` + `Bonds`).
- `--direction {x,y,z,0,1,2}` (default `x`): the loading axis OBC weights against.
- `--validate-stock`: also runs stock `igraph.Graph.edge_betweenness()` and fails if
  native GEBC disagrees (`rtol=1e-12, atol=1e-8`) - good for a first run on new data,
  skip it afterward to avoid a second full traversal.
- `--inspect`: parses and validates the input, prints its topology/box as JSON, and
  stops - no centrality is computed.
- `--sample-vertices N`: fast debug run on an induced neighborhood of `N` vertices
  (a genuinely different, smaller graph - not an estimate of the full network).
- `--dbc-length L`: raw DBC is always distance-only (no automatic length division);
  this selects a divisor for the *postprocessed* `DBC` output only. Without it, `DBC`
  carries only pair normalization. Apply any length later via
  `CentralityResult.dbc_for_length(L)` without recomputing.
- `--output-file` (default `<input>.<axis>.spatial.npz`) and a matching `.json`
  report; `--overwrite` to replace existing output, `--quiet` to suppress progress.

Output arrays (`bond_ids`, `edges`, `raw_GEBC_igraph`, `raw_OBC`, `raw_DBC`, `OBC`,
`DBC`, ...) are index-aligned with the input's bond rows - see
`../vendor/python-igraph/analysis/README.md` for the full output/scientific-convention
reference and the extensive validation already recorded there.

## As a library

```python
from graph_analysis import PolymerNetwork, NetworkAnalysisPipeline, PipelineConfig

network = PolymerNetwork.from_lammps("atrp_s3_sys3_0.data")  # a full LAMMPS data file
# network = PolymerNetwork.from_bond_list("edges_0.bds")     # legacy: a pre-trimmed bond table
pipeline = NetworkAnalysisPipeline(network, PipelineConfig(sys_name="run1", ncpus=4))
result = pipeline.run()  # a ScanlanCaseResult: result.refined, result.raw

# Needs the in-house igraph build (../vendor/) - run with vendor/.spatial-venv/bin/python
from graph_analysis.betweenness import SpatialBetweennessAnalysis, AnalysisConfig

analysis = SpatialBetweennessAnalysis.from_lammps("atrp_s1_sys1_0.data", AnalysisConfig(direction="x"))
result = analysis.run()  # a CentralityResult: result.raw_gebc, result.obc, result.dbc, ...
```

See `graph_analysis/__init__.py` for the topology/shortest-path public API
(`PolymerNetwork`, `CrosslinkAction`, `CrosslinkPathAnalyzer`, `LoopAnalyzer`,
`LoopBasis`, `NetworkStats`, `RefinedStats`, `ScanlanCaseAnalyzer`,
`ScanlanCaseResult`, `NetworkAnalysisPipeline`, `PipelineConfig`) and
`graph_analysis/betweenness/__init__.py` for the betweenness public API
(`SpatialBetweennessAnalysis`, `AnalysisConfig`, `Direction`, `PolymerNetwork`,
`LammpsDataReader`, `CentralityResult`, `ValidationReport`) - deliberately not
re-exported at the top level since its `PolymerNetwork` is a distinct, array-based
representation, not interchangeable with the topology one.
