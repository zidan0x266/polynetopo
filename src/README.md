# Polymer network topology analysis

A single CLI (`graph_analysis.cli`, exposed here as `main.py`) with subcommands for
the crosslinked-network analysis pipeline: elastically effective strands, dangling
chains, primary/higher-order loops, and the Scanlan-Case crosslinker density
refinement (used for the ATRP/FRP network papers in `../references`).

Run `python3 main.py --help` (or `--help` on any subcommand) for the full option list.

## One-shot (recommended for a new project)

`analyze` runs every stage - strands, danglings, loops, Scanlan-Case refinement -
from a single topology load, with no intermediate files:

```bash
python3 main.py analyze --input-file ../examples/atrp/edges_0.bds \
    --output-file results.pkl --raw --ncpus 4 --max-chunk-size 10000
```

Writes `results.pkl` (a `graph_analysis.RefinedStats`: the Scanlan-Case-refined
strands/danglings/loops plus the crosslinker density it produced). Add `--raw` to
also write `results_raw.pkl` (a `graph_analysis.NetworkStats`: the pre-refinement
strands/danglings/loops, for sanity-checking a specific configuration without
redoing the analysis). `--nvertices` is optional - if omitted it's inferred from the
topology (see the reproducibility note in `graph_analysis/loops.py` before relying on
this to reproduce an old run bit-for-bit: vertex count doesn't affect the result, but
edge *order* does, which is why `--topology rmc`/`lmp` both preserve the file's
original line order rather than going through `networkx`'s own edge iteration).

## Step-by-step (the original workflow, still available)

Useful for inspecting an intermediate stage, or resuming after a crash - each stage
caches its crosslinker search matrix under `cache-<sys-name>-*/`.

```bash
# Strands: shortest paths between crosslinkers (degree >= 3)
python3 main.py distances --sys-name strands --input-file ../examples/atrp/edges_0.bds \
    --output-file strands.pkl --action crosslinkers --max-chunk-size 10000 --ncpus 4 --topology lmp
rm -r cache-*

# Dangling chains: shortest paths from a crosslinker to a chain end (degree == 1)
python3 main.py distances --sys-name danglings --input-file ../examples/atrp/edges_0.bds \
    --output-file danglings.pkl --action danglings --ncpus 4 --max-chunk-size 10000 --topology lmp
rm -r cache-*

# Minimum cycle basis (loops), with per-loop crosslinker counts
python3 main.py loops --input-file ../examples/atrp/edges_0.bds --output-file igloops.pkl --nvertices 106250

# Scanlan-Case refinement of the strands/danglings/loops above (also takes --raw)
python3 main.py scanlan-case --input-file ../examples/atrp/edges_0.bds \
    --strands-file strands.pkl --danglings-file danglings.pkl --loops-file igloops.pkl \
    --output-file Scanlan_Case.pkl

# Compare the entry counts of two pickled result files
python3 main.py compare strands.pkl danglings.pkl
```

Outputs land in `../examples/atrp/*.pkl`.

`Scanlan_Case.pkl`/`results.pkl` now pickle a `graph_analysis.RefinedStats` dataclass
rather than the old positional nested-list structure; call `.to_legacy()` on the
`ScanlanCaseAnalyzer.run()` return value (before pickling) to recover the previous
`[[refined...], [raw...], [counts...]]` shape if needed for older downstream scripts.

## As a library

```python
from graph_analysis import PolymerNetwork, NetworkAnalysisPipeline, PipelineConfig

network = PolymerNetwork.from_lammps("edges_0.bds")
pipeline = NetworkAnalysisPipeline(network, PipelineConfig(sys_name="run1", ncpus=4))
result = pipeline.run()  # a ScanlanCaseResult: result.refined, result.raw
```

See `graph_analysis/__init__.py` for the full public API (`PolymerNetwork`,
`CrosslinkAction`, `CrosslinkPathAnalyzer`, `LoopAnalyzer`, `LoopBasis`, `NetworkStats`,
`RefinedStats`, `ScanlanCaseAnalyzer`, `ScanlanCaseResult`, `NetworkAnalysisPipeline`,
`PipelineConfig`).
