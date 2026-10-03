# PolyNetopo

Topology, shortest-path and spatial betweenness analysis for crosslinked polymer
networks: elastically effective strands, dangling chains, primary/higher-order loops,
Scanlan-Case crosslinker density refinement, and native GEBC/OBC/DBC edge betweenness.
Built on `networkx` (shortest-path search), stock `igraph` (minimum cycle basis), and
an in-house `igraph` build (spatial betweenness), with a multiprocessing-accelerated
search over large networks.

This code produced the analysis in the three papers cited at the bottom of this README.

## What it does

Given a polymer network topology - a full LAMMPS data file by default, or a legacy
pre-trimmed bond table / RMC-style edge list - `--function topology` identifies:

- **Elastically effective strands** - the shortest simple paths between crosslinkers
  (nodes of degree >= 3).
- **Dangling chains** - shortest simple paths from a crosslinker to a free chain end
  (degree == 1); topological defects that don't bear load.
- **Loops** - the network's minimum cycle basis, with a crosslinker count per loop.
  A loop with exactly one crosslinker on it is a *primary loop*.
- **Scanlan-Case-refined crosslinker density** - reclassifies crosslinkers that carry a
  primary loop or dangling chain(s) as linear units rather than real network junctions,
  folding their strands/danglings together and recomputing effective vs. total
  crosslinker counts.

`--function betweenness` computes, for every bond in one pass:

- **GEBC** - geodesic edge betweenness centrality, the standard topological metric.
- **OBC** - orientation-weighted betweenness, by alignment with a loading axis.
- **DBC** - distance-weighted betweenness, by source-target spatial separation.

## Installation

```bash
pip install -e .
```

Requires Python >= 3.10, `networkx`, `numpy`, `igraph` (the `python-igraph` package),
and `tqdm`. `--function betweenness` additionally needs the in-house igraph build in
`vendor/` - see [`vendor/README.md`](vendor/README.md).

## Quick start

Run the whole topology pipeline - strands, danglings, loops, and their Scanlan-Case
refinement - from a single topology load, with no intermediate files:

```bash
python3 src/main.py --function topology \
    --input-file examples/shortest_path/atrp_s3_sys3_0.data \
    --output-file results.pkl --raw --ncpus 4 --max-chunk-size 10000
```

`--input-file` here is a full LAMMPS data file (`write_data` output) - the default,
primary format (`--format lammps`); only its Bonds section is read. A pre-trimmed bare
bond table is also supported as a legacy format (`--format bond-list`), e.g.
`examples/atrp/edges_0.bds`.

This writes `results.pkl` (a `graph_analysis.RefinedStats`: the refined
strands/danglings/loops plus the crosslinker density they produced) and, since `--raw`
was passed, `results_raw.pkl` (a `graph_analysis.NetworkStats`: the pre-refinement data,
for sanity-checking a configuration without redoing the analysis). `--stage
loops`/`scanlan-case` runs just one stage - useful to resume after a crash, or to run a
memory-heavy stage on a bigger machine separately.

See [Using `--function betweenness` correctly](#using---function-betweenness-correctly)
below for spatial betweenness.

As a library:

```python
from graph_analysis import PolymerNetwork, NetworkAnalysisPipeline, PipelineConfig

network = PolymerNetwork.from_lammps("atrp_s3_sys3_0.data")  # a full LAMMPS data file
pipeline = NetworkAnalysisPipeline(network, PipelineConfig(sys_name="run1", ncpus=4))
result = pipeline.run()  # a ScanlanCaseResult: result.refined, result.raw

from graph_analysis.betweenness import SpatialBetweennessAnalysis  # needs the in-house build
result = SpatialBetweennessAnalysis.from_lammps("examples/betweenness/atrp_s1_sys1_0.data").run()
```

See [`src/README.md`](src/README.md) for the full CLI reference and complete public API.

## Package layout

```
src/graph_analysis/
    network.py        PolymerNetwork - loads a topology (LAMMPS data file, or legacy formats)
    crosslink.py       CrosslinkPathAnalyzer - shortest simple paths between crosslinkers
    loops.py           LoopAnalyzer - minimum cycle basis loop detection
    scanlan_case.py    ScanlanCaseAnalyzer - crosslinker density refinement
    pipeline.py        NetworkAnalysisPipeline - runs all of the above in one pass
    betweenness/       SpatialBetweennessAnalysis - native GEBC/OBC/DBC (needs vendor/)
    cli.py             command-line entry point (--function topology/shortest-path/betweenness)
vendor/
    igraph/            git submodule: in-house igraph C core (spatial-betweenness-1.0.0)
    python-igraph/     git submodule: in-house Python bindings (same branch)
```

Every graph view used by `--function topology`/`shortest-path` (the `networkx.Graph`
used for strand/dangling search and degree lookups, and the `igraph.Graph` used for
loop detection) is built from the same single, file-order-preserving parse in
`PolymerNetwork` - so node IDs (raw atom IDs from the topology file) are guaranteed
consistent across strands, danglings, and loops, and a loop's member atoms can be
cross-referenced against strand/dangling endpoints directly by ID.

## Using `--function betweenness` correctly

The three most common mistakes are using the wrong Python interpreter, skipping
validation on new data, and misreading the raw DBC convention - this section covers
all three.

### 1. Build once per machine

```bash
git submodule update --init vendor/igraph vendor/python-igraph
echo "1.0.0" > vendor/igraph/IGRAPH_VERSION   # the fork has no git tags; see vendor/README.md
cd vendor && bash python-igraph/scripts/build_spatial.sh && cd ..
vendor/.spatial-venv/bin/python -m pip install --no-build-isolation --no-deps networkx tqdm
vendor/.spatial-venv/bin/python -m pip install -e .
```

Takes ~15-20 minutes (full C+Python build and test suite); do this once per machine -
binaries from one machine/OS don't work on another. Full details, including the exact
build error you'd hit without the `IGRAPH_VERSION` workaround, are in
[`vendor/README.md`](vendor/README.md).

### 2. Always use that build's Python, never your regular one

```bash
vendor/.spatial-venv/bin/python -c 'import igraph as ig; print(hasattr(ig.Graph, "edge_betweenness_spatial"))'
# must print True - if it prints False or raises ImportError, you're using the wrong interpreter
```

Every command below uses `vendor/.spatial-venv/bin/python`, not plain `python3`. Stock
PyPI `igraph` (what `--function topology`/`shortest-path` use) doesn't have
`edge_betweenness_spatial`; calling `--function betweenness` with it raises a clear
`RuntimeError` rather than silently falling back to something else.

### 3. Inspect new input before computing anything

```bash
vendor/.spatial-venv/bin/python src/main.py --function betweenness \
    --input-file examples/betweenness/atrp_s1_sys1_0.data --inspect
```

Parses and validates the file, prints its topology/box as JSON, and does no
centrality calculation. Confirms atom/bond counts and the *actual* box lengths before
you commit to a multi-minute run - a deformed network's box is not the cubic
undeformed one, and OBC/DBC depend on it.

### 4. Validate on first use with new data, then drop `--validate-stock` for routine runs

```bash
# First run on a new network: also runs stock GEBC and fails loudly if native GEBC disagrees
vendor/.spatial-venv/bin/python src/main.py --function betweenness \
    --input-file examples/betweenness/atrp_s1_sys1_0.data --direction x --validate-stock
```

`--validate-stock` roughly doubles runtime (a second full traversal) - worth it once
per new network/machine combination, skip it for routine repeated runs afterward.
For a quick sanity check before a full run, use `--sample-vertices 256
--validate-stock` first (seconds, not minutes - but a genuinely different, smaller
graph, not an estimate of the full network).

### 5. Read the output correctly - raw DBC has no automatic length division

Output is a pickled `CentralityResult` - a plain pickle of a dataclass, like every
other output in this package, not a NumPy archive:

```python
import pickle
with open("examples/betweenness/atrp_s1_sys1_0.x.spatial.pkl", "rb") as f:
    result = pickle.load(f)  # needs graph_analysis importable, e.g. pip install -e .
raw_gebc, raw_obc, raw_dbc = result.raw_gebc, result.raw_obc, result.raw_dbc
obc, dbc = result.obc, result.dbc  # pair-normalized properties; dbc also length-divided IF --dbc-length was passed
```

`raw_dbc` is **always** the distance-only sum (no `Lx`/bond-length denominator) -
apply any reference length yourself, as many times as you like, without recomputing:
`result.dbc_for_length(L)` or `result.dbc_for_length(L, normalized=True)`. Passing
`--dbc-length L` at the CLI only changes the *postprocessed* `dbc` property; `raw_dbc`
is unaffected either way. See `vendor/python-igraph/analysis/README.md` for the full
convention reference (normalization factors, what each field indexes, migration notes
for files saved before this convention) - written against that package's original
NumPy-archive field names, which map directly onto this pickle's attributes (e.g.
`raw_GEBC_igraph`/`OBC`/`DBC` there are `result.raw_gebc`/`obc`/`dbc` here).

## Example and test datasets

Each folder has its own `TUTORIAL.md` - a step-by-step walkthrough of the
`--function` it's named after, using that folder's data as the worked example.

- [`examples/atrp/`](examples/atrp/TUTORIAL.md) (`--function topology`) - an
  ATRP-synthesized network (`edges_0.bds`, a legacy pre-trimmed bond table, 103,149
  atoms after removing unreacted monomer at 96% conversion) and its topology analysis
  outputs (`strands.pkl`, `danglings.pkl`, `igloops.pkl`, `Scanlan_Case.pkl`). Also
  the tutorial for the legacy bond-list format.
- [`examples/shortest_path/`](examples/shortest_path/TUTORIAL.md)
  (`--function shortest-path`) - a denser network (`atrp_s3_sys3_0.data`, a full
  LAMMPS data file, 104,704 atoms, 9,537 crosslinkers) with the same set of topology
  outputs. Also the tutorial for the default full-data-file format.
- [`examples/betweenness/`](examples/betweenness/TUTORIAL.md)
  (`--function betweenness`) - the undeformed (`atrp_s1_sys1_0.data`) and a deformed
  (`atrp_s1_sys1_42.data`) state of a 103,149-atom/105,869-bond network, used to
  validate `--function betweenness` against stock GEBC and independent path
  enumeration (see `vendor/python-igraph/SPATIAL_BETWEENNESS.md`).

All reproducible with the commands in `src/README.md` or each folder's tutorial.

## Reproducibility notes

- **Node IDs are raw atom IDs everywhere** - the topology file's own numbering, used
  consistently as `networkx` node labels, `igraph` vertex indices, and every strand/
  dangling/loop entry. No remapping step is needed to cross-reference between them.
- **Loop identity depends on edge order, not vertex count.** A network's minimum cycle
  basis isn't always unique; when equally-short candidate cycles tie, `igraph`'s
  tie-breaking depends on the order edges were added, not on how many unused vertex IDs
  the graph was padded out to (`--nvertices`, if given, is otherwise safe to omit or
  over-specify). See `src/graph_analysis/loops.py` for the measurements behind this.
- **A latent bug in the original strand-merging code was found and fixed** during this
  package's development: `analysis_Scanlan_Case.py`'s merge loop removed items from a
  list while iterating over it directly, silently skipping the element after any
  removed one. On this repo's example dataset that left 13 of 4632 refined strands
  (~0.3%) merged incompletely; aggregate crosslinker-density statistics are unaffected.
  See `src/graph_analysis/scanlan_case.py` for details.

## Citing this work

If you use this code, please cite:

> Zhang, Z.; Krajniak, J.; Dookhith, A. Z.; Tian, Y.; Sachar, H. S.; Marioni, N.;
> Duncan, T. J.; Liu, J.; Sanoja, G. E.; Ganesan, V. Topology and Mechanical Properties
> of Polymer Networks Formed under Free Radical and Atom Transfer Radical
> Polymerizations. *Macromolecules* **2025**.
> DOI: [10.1021/acs.macromol.4c03138](https://doi.org/10.1021/acs.macromol.4c03138)

> Zhang, Z.; Krajniak, J.; Dookhith, A. Z.; Zhang, H.; Tian, Y.; Sachar, H. S.; Marioni,
> N.; Duncan, T. J.; Liu, J.; Sanoja, G. E.; Ganesan, V. Fracture Behavior of
> Monodisperse and Polydisperse Polymer Networks: Effects of Topology and
> Strand-Length Dispersity. *Macromolecules* **2025**.
> DOI: [10.1021/acs.macromol.5c01715](https://doi.org/10.1021/acs.macromol.5c01715)

> Bai, X.; Ding, Z.; Qi, S.; Peng, L.; Huang, X.; Yan, D.; Man, X. Spatially Weighted
> Betweenness Centrality for Predicting Fracture in End-Linked Polymer Networks.
> *Macromolecules* **2026**, *59*, 10294-10303.
> DOI: [10.1021/acs.macromol.6c00974](https://doi.org/10.1021/acs.macromol.6c00974)

```bibtex
@article{zhang2025topology,
  title   = {Topology and Mechanical Properties of Polymer Networks Formed under Free Radical and Atom Transfer Radical Polymerizations},
  author  = {Zhang, Zidan and Krajniak, Jakub and Dookhith, Aaliyah Z. and Tian, Yuan and Sachar, Harnoor S. and Marioni, Nico and Duncan, Tyler J. and Liu, Jun and Sanoja, Gabriel E. and Ganesan, Venkat},
  journal = {Macromolecules},
  year    = {2025},
  doi     = {10.1021/acs.macromol.4c03138}
}

@article{zhang2025fracture,
  title   = {Fracture Behavior of Monodisperse and Polydisperse Polymer Networks: Effects of Topology and Strand-Length Dispersity},
  author  = {Zhang, Zidan and Krajniak, Jakub and Dookhith, Aaliyah Z. and Zhang, Han and Tian, Yuan and Sachar, Harnoor S. and Marioni, Nico and Duncan, Tyler J. and Liu, Jun and Sanoja, Gabriel E. and Ganesan, Venkat},
  journal = {Macromolecules},
  year    = {2025},
  doi     = {10.1021/acs.macromol.5c01715}
}

@article{bai2026gebc,
  title   = {Spatially Weighted Betweenness Centrality for Predicting Fracture in End-Linked Polymer Networks},
  author  = {Bai, Xianjue and Ding, Zhaoyu and Qi, Shuanhu and Peng, Li and Huang, Xianbo and Yan, Dadong and Man, Xingkun},
  journal = {Macromolecules},
  year    = {2026},
  volume  = {59},
  pages   = {10294--10303},
  doi     = {10.1021/acs.macromol.6c00974}
}
```
