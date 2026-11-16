# Polymer network topology, shortest-path and betweenness analysis

A single CLI (`graph_analysis.cli`, exposed here as `main.py`) covering the
crosslinked-network analysis pipeline: elastically effective strands, dangling
chains, primary/higher-order loops, Scanlan-Case crosslinker density refinement
(used for the ATRP/FRP network papers), shortest paths across the periodic boundary
(Yu & Jackson, which govern bond rupture), and native GEBC/OBC/DBC spatial edge
betweenness (used for the fracture-prediction paper).

One flag, `--function {topology, shortest-path, betweenness, compare}` (default
`topology`), selects what to run; run `python3 main.py --help` for the full option
list, grouped by which function each belongs to. The commands below use the shared
example network `../examples/data/atrp_s1_sys1_0.data`.

## Options at a glance

| Option | Functions | Meaning |
|---|---|---|
| `--function` | all | `topology` (default), `shortest-path`, `betweenness`, `compare` |
| `--input-file` | topology, shortest-path, betweenness | the network, normally a full LAMMPS data file |
| `--output-file` | all | output pickle; each function's default is given below |
| `--format` | topology | `lammps` (default), `bond-list`, `rmc` (see below) |
| `--direction {x,y,z,0,1,2}` | betweenness, shortest-path | betweenness: the loading axis OBC weights against; shortest-path: the periodic axis the paths wrap along (default `x`) |
| `--overwrite` | betweenness, shortest-path | replace an existing output file |
| `--stage` | topology | `full` (default), `strands`, `danglings`, `loops`, `scanlan-case` |
| `--raw` | topology | also write the pre-refinement result (`full`, `scanlan-case`) |
| `--nvertices N` | topology | number of nodes; inferred if omitted |
| `--strands-file`, `--danglings-file`, `--loops-file` | topology `scanlan-case` | the three files it refines |
| `--max-period M` | shortest-path | search periods p = 1..M (default 8, the paper's M) |
| `--buffer-images N` | shortest-path | extra images on each side (default 1) |
| `--slice-max L` | shortest-path | only start paths from atoms within L of the box's low face |
| `--keep-redundant` | shortest-path | skip the removal of redundant paths (Cases 1 and 2) |
| `--min-sp-only` | shortest-path | also apply the paper's optional Case 3 (only the shortest SP of each atom) |
| `--ncpus N` | shortest-path | worker processes for the per-atom search (default 1) |
| `--checkpoint-file F` | shortest-path | save progress after each period, resume if it exists |
| `--validate-stock`, `--inspect`, `--sample-vertices N`, `--dbc-length L`, `--quiet` | betweenness | see below |
| `--file1`, `--file2` | compare | print the entry counts of two pickled results |

## Topology formats

`--format` (default `lammps`) selects how `--input-file` is read for `--function topology`:

- `lammps` (default) - a full LAMMPS data file (`write_data` output, with its
  Masses/Atoms/Velocities/Bonds/... sections); only the Bonds section is read.
- `bond-list` - **legacy**: a bare bond table, `bond-id bond-type atom_i atom_j` per
  line.
- `rmc` - **legacy**: an RMC-style whitespace edge list.

`--function shortest-path` and `--function betweenness` always read a full LAMMPS data
file directly (atoms with coordinates, bonds, and an orthorhombic box; see
`graph_analysis/lammps.py`).

## `--function topology`

`--stage full` (the default) runs every stage - strands, danglings, loops,
Scanlan-Case refinement - from a single topology load:

```bash
python3 main.py --function topology \
    --input-file ../examples/data/atrp_s1_sys1_0.data \
    --output-file results.pkl --raw
```

It writes, next to each other:

- `results.pkl` - a `graph_analysis.RefinedStats`: the Scanlan-Case-refined
  strands/danglings/loops plus the crosslinker density they produced.
- `results_strands.pkl`, `results_danglings.pkl` - flat lists of atom-ID paths.
- `results_loops.pkl` - `[loops, crosslinker_counts]`.
- with `--raw`, `results_raw.pkl` - a `graph_analysis.NetworkStats`: the pre-refinement
  strands/danglings/loops in one object.

`--nvertices` is optional - if omitted it's inferred from the topology.

### Single stages

Each stage writes only its own output and never triggers the refinement, because the
three files it needs may not exist. Useful for inspecting an intermediate stage, running
the memory-heavy loop search on a bigger machine separately, or rerunning the
refinement:

```bash
D=../examples/data/atrp_s1_sys1_0.data
python3 main.py --function topology --stage strands   --input-file $D --output-file strands.pkl
python3 main.py --function topology --stage danglings --input-file $D --output-file danglings.pkl
python3 main.py --function topology --stage loops     --input-file $D --output-file igloops.pkl

# Scanlan-Case refinement of the three files above (also takes --raw); the files written
# by `--stage full` (results_strands.pkl, ...) work the same way
python3 main.py --function topology --stage scanlan-case --input-file $D \
    --strands-file strands.pkl --danglings-file danglings.pkl --loops-file igloops.pkl \
    --output-file Scanlan_Case.pkl

# Compare the entry counts of two pickled result files
python3 main.py --function compare --file1 strands.pkl --file2 danglings.pkl
```

`--stage loops` takes only `--input-file`, `--output-file`, and the optional `--format`
and `--nvertices`; it has no options controlling the search itself (it calls igraph's
`minimum_cycle_basis()`). In a full run on the 103,149-atom example, on a laptop, the
loop search took about 2.4 minutes and the whole process peaked near 8.6 GB.

`Scanlan_Case.pkl`/`results.pkl` pickle a `graph_analysis.RefinedStats` dataclass; call
`.to_legacy()` on the `ScanlanCaseAnalyzer.run()` return value (before pickling) for the
nested-list form `[[refined...], [raw...], [counts...]]`.

## `--function shortest-path`

The method of Yu & Jackson, "Shortest Paths Govern Bond Rupture in Thermoset Networks"
(Macromolecules 2025), from the authors' code. For each period p = 1..M it builds a
graph of p + 1 periodic images (plus `--buffer-images` on each side), finds for every
atom the shortest path from itself in image 0 to itself in image p, and records its
length divided by p: **D^g = edges / p**. The README's "What it does" section gives the
definitions of valid paths and of Cases 1, 2 and 3.

```bash
python3 main.py --function shortest-path \
    --input-file ../examples/data/atrp_s1_sys1_0.data --direction x \
    --max-period 4 --ncpus 8
```

- `--max-period M` (default 8): periods 1..M inclusive. The search stops earlier at the
  first period without a valid path.
- `--buffer-images N` (default 1): extra images on each side, so a path can leave the
  stack between its ends; `0` keeps every path inside it.
- `--slice-max L`: start paths only from atoms within L of the box's low face along
  `--direction`. This reduces cost; a quick trial is `--slice-max 1.0 --max-period 2`.
- `--keep-redundant`: skip the removal of redundant paths (Cases 1 and 2).
- `--min-sp-only`: also apply Case 3, keeping only the shortest SP of each start atom
  (ties go to the first; the result is ordered by start atom). Off by default.
- `--ncpus N`: worker processes, each holding one copy of the multi-image graph.
- `--checkpoint-file F`: write progress after each period; rerunning the same command
  resumes after the last finished period. `--max-period` and the redundancy options may
  differ between runs, the network, direction, `--buffer-images` and `--slice-max` may
  not.
- `--output-file` (default `<input>.<axis>.sp.pkl`); `--overwrite` to replace it.

The input may not contain an atom pair bonded twice or an atom bonded to itself; this
raises an error naming the pair.

The output is a pickled `graph_analysis.shortest_path.ShortestPathResult`:

```python
import pickle
with open("atrp_s1_sys1_0.x.sp.pkl", "rb") as f:
    sp = pickle.load(f)               # needs graph_analysis importable, e.g. PYTHONPATH=src
sp.lengths                            # D^g of every path, aligned with sp.paths
sp.paths[0]                           # (n, 2) rows of (atom ID, image); first and last atom are equal
sp.min_length                         # D^g_min, the first (shortest) SP
sp.min_length_by_period               # rows of (p, smallest D^g found at p), before Cases 1-3
sp.summary()                          # the dictionary the command prints
```

Images are counted from image 0, the box as given: a path starts at image 0 and ends
at image p, and an image of -1, or above p, is a buffer image it passed through. Atom
IDs are the raw IDs of the data file. `min_length_by_period` is computed before
redundant paths are removed, so it is the same with or without `--min-sp-only`.

## `--function betweenness`

Native GEBC/OBC/DBC edge betweenness, computed together in one pass. Needs the
in-house igraph build in `../vendor/` (see `../vendor/README.md`) - run with that
build's own Python, not your regular environment's:

```bash
../vendor/.spatial-venv/bin/python main.py --function betweenness \
    --input-file ../examples/data/atrp_s1_sys1_0.data --direction x --validate-stock
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
- `--output-file` (default `<input>.<axis>.spatial.pkl`); `--overwrite` to replace
  existing output, `--quiet` to suppress progress.

The output is a pickled `graph_analysis.betweenness.CentralityResult`:

```python
import pickle
with open("atrp_s1_sys1_0.x.spatial.pkl", "rb") as f:
    result = pickle.load(f)  # needs graph_analysis importable, e.g. PYTHONPATH=src
bond_ids = result.network.bond_ids       # aligned with every array below
raw_gebc, raw_obc, raw_dbc = result.raw_gebc, result.raw_obc, result.raw_dbc
obc, dbc = result.obc, result.dbc        # pair-normalized properties, not stored separately
```

See `../vendor/python-igraph/analysis/README.md` for the full convention reference.

## As a library

```python
from graph_analysis import PolymerNetwork, NetworkAnalysisPipeline

network = PolymerNetwork.from_lammps("atrp_s1_sys1_0.data")  # a full LAMMPS data file
# network = PolymerNetwork.from_bond_list("edges_0.bds")     # legacy: a pre-trimmed bond table
pipeline = NetworkAnalysisPipeline(network)
result = pipeline.run()  # a ScanlanCaseResult: result.refined, result.raw

from graph_analysis.shortest_path import PeriodicShortestPathAnalysis, ShortestPathConfig

config = ShortestPathConfig(direction="x", max_period=4, ncpus=8)  # min_sp_only=True for Case 3
sp = PeriodicShortestPathAnalysis.from_lammps("atrp_s1_sys1_0.data", config).run()

# Needs the in-house igraph build (../vendor/) - run with vendor/.spatial-venv/bin/python
from graph_analysis.betweenness import SpatialBetweennessAnalysis, AnalysisConfig

analysis = SpatialBetweennessAnalysis.from_lammps("atrp_s1_sys1_0.data", AnalysisConfig(direction="x"))
result = analysis.run()  # a CentralityResult: result.raw_gebc, result.obc, result.dbc, ...
```

Public API by module:

- `graph_analysis` (topology): `PolymerNetwork`, `CrosslinkAction`,
  `CrosslinkPathAnalyzer`, `LoopAnalyzer`, `LoopBasis`, `NetworkStats`, `RefinedStats`,
  `ScanlanCaseAnalyzer`, `ScanlanCaseResult`, `NetworkAnalysisPipeline`.
- `graph_analysis.lammps` (shared by the spatial analyses): `LammpsDataReader`,
  `LammpsFormatError`, `LammpsNetwork` (atoms, coordinates, bonds and box), `Direction`.
- `graph_analysis.shortest_path`: `PeriodicShortestPathAnalysis`, `ShortestPathConfig`,
  `ShortestPathResult`.
- `graph_analysis.betweenness`: `SpatialBetweennessAnalysis`, `AnalysisConfig`,
  `CentralityResult`, `ValidationReport`, `CentralityValidationError`, and re-exports of
  the `lammps` names above.
