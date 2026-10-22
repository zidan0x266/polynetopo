# PolyNetopo

Topology, shortest-path and spatial betweenness analysis for crosslinked polymer
networks, from a LAMMPS data file:

| `--function` | What it computes | Needs |
|---|---|---|
| `topology` | elastically effective strands, dangling chains, primary/higher-order loops, and the Scanlan-Case crosslinker density refinement | bonds only |
| `shortest-path` | the shortest paths that wrap through the periodic box, which govern where and when the first bond breaks (Yu & Jackson, 2025) | bonds, coordinates, box |
| `betweenness` | native GEBC / OBC / DBC edge betweenness | bonds, coordinates, box |

This code produced the analysis in the first three papers cited at the bottom of this
README. Its shortest-path analysis implements the method of the fourth, Yu & Jackson,
ported from the authors' code
([TheJacksonLab/shortest_path_thermoset](https://github.com/TheJacksonLab/shortest_path_thermoset)).
It is built on `networkx` and `igraph`; the betweenness analysis needs the in-house
`igraph` build in `vendor/`, which also runs every other function.

## What it does

### `--function topology`

Given a network topology - a full LAMMPS data file by default, or a legacy bare bond
table / RMC-style edge list (`--format bond-list` / `rmc`) - it identifies:

- **Elastically effective strands** - the chains of degree-2 atoms between two
  crosslinkers (nodes of degree >= 3). Every distinct chain is found, including several
  parallel chains between the same two crosslinkers.
- **Dangling chains** - chains from a crosslinker to a free chain end (degree == 1);
  topological defects that don't bear load.
- **Loops** - the network's minimum cycle basis, with a crosslinker count per loop.
  A loop with exactly one crosslinker on it is a *primary loop*.
- **Scanlan-Case-refined crosslinker density** - reclassifies crosslinkers that carry a
  primary loop or dangling chain(s) as linear units rather than real network junctions,
  folding their strands/danglings together and recomputing effective vs. total
  crosslinker counts.

`--stage` picks what runs:

| `--stage` | Runs | Writes |
|---|---|---|
| `full` (default) | strands, danglings, loops, then Scanlan-Case, from one load | `results.pkl` (refined result), `results_strands.pkl`, `results_danglings.pkl`, `results_loops.pkl`; with `--raw` also `results_raw.pkl` |
| `strands` | only the strand search | a flat list of atom-ID paths |
| `danglings` | only the dangling search | a flat list of atom-ID paths |
| `loops` | only the minimum cycle basis | `[loops, crosslinker_counts]` |
| `scanlan-case` | only the refinement, from `--strands-file`, `--danglings-file` and `--loops-file` | a refined result |

A single stage never runs the refinement, since the three files it needs may not exist.
`full` writes them next to the refined result, so the refinement can be rerun without
redoing the slow loop search:

```bash
python3 src/main.py --function topology --stage scanlan-case \
    --input-file examples/data/atrp_s1_sys1_0.data \
    --strands-file results_strands.pkl --danglings-file results_danglings.pkl \
    --loops-file results_loops.pkl --output-file results_again.pkl
```

### `--function shortest-path`

The method of Yu & Jackson, "Shortest Paths Govern Bond Rupture in Thermoset Networks".
A shortest path (SP) that wraps through the periodic box can be longer than one image,
so for each period p = 1, 2, ..., M (`--max-period`, the paper's M, default 8):

1. **Multi-image graph.** p + 1 periodic images of the network along `--direction`
   (image 0 is the box as given), plus `--buffer-images` extra images on each side
   (default 1). Bonds stay inside an image, and bonds that cross the periodic boundary
   connect neighbouring images. The buffer images let a path leave the stack between
   its ends instead of being forced to stay inside it; 0 keeps them inside.
2. **Shortest paths.** For every atom, the shortest path from that atom in image 0 to
   the same atom in image p. Its length is normalized by p, **D^g = edges / p**, the
   number of bonds per box length, so paths that wrap different numbers of times are
   comparable. A path is *valid* only if the start/end atom is the only atom it visits
   twice; any other path is discarded (it only repeats shorter ones).
3. **Redundant paths** are removed when an already kept path of equal D^g contains all
   of a path's atoms. This covers the paper's Case 1 (the same path found from
   different start atoms) and, together with the validity rule, Case 2 (a path at larger
   p that repeats one at smaller p). `--keep-redundant` skips it.

The search stops at the first period with no valid path. The smallest D^g is the first
(shortest) SP; the paper uses it to predict where the first bond breaks and, through
D^g_min / L0, the strain at which it breaks.

Optional **Case 3** (`--min-sp-only`): of the SPs starting at one atom, keep only the
shortest. It is off by default, so all valid paths of an atom are kept.

```bash
python3 src/main.py --function shortest-path \
    --input-file examples/data/atrp_s1_sys1_0.data --direction x --ncpus 8 \
    --output-file sp.pkl
```

The output is a pickled `ShortestPathResult`: `lengths` (D^g per path), `paths` (one
`(n, 2)` array of `(atom ID, image)` rows per path, images counted from image 0, so
`-1` or above p is a buffer image), `min_length_by_period`, and `summary()`. The input
must have no repeated bonds (an atom pair bonded twice) - a path lists atoms, so it
could not say which of two parallel bonds it used; this raises an error naming the
pair.

Cost: one breadth-first search per start atom over a graph of p + 1 + 2 x
`--buffer-images` copies of the network, so time grows with atoms x images x
(atoms + bonds), and memory with p (a path has about D^g x p atoms). `--ncpus` splits the
start atoms over worker processes (each holds one copy of the graph). To try a network
cheaply, restrict the start atoms with `--slice-max 1.0` and use a small `--max-period`;
`--checkpoint-file` saves progress after each period so a long run can resume.

### `--function betweenness`

Computes, for every bond in one pass:

- **GEBC** - geodesic edge betweenness centrality, the standard topological metric.
- **OBC** - orientation-weighted betweenness, by alignment with a loading axis.
- **DBC** - distance-weighted betweenness, by source-target spatial separation.

See [Using `--function betweenness` correctly](#using---function-betweenness-correctly)
below.

## Note on strand counts from earlier releases

Earlier releases found strands with a pairwise Yen's-algorithm search that missed
some strands when three or more separate chains join the same two crosslinkers (its
`ignore_nodes` set was shared across iterations, so the source node blocked later
searches). The search is now a chain walk that finds every chain. On the 104,704-atom
sys3 network this adds 57 strands (17,934 -> 17,991, ~0.3%) and on the 103,149-atom sys1
network 18 (5,138 -> 5,156); dangling chains and loops are unchanged. Strand-derived
numbers from earlier runs, including the Scanlan-Case refinement, may differ slightly.

## Installation

Requires Python >= 3.10 with `networkx` (>= 3.0), `numpy`, `igraph` (the `python-igraph`
package) and `tqdm`. Install with `pip install -e .`, or skip the install and run
`python3 src/main.py` from the repository root, which works without it.
`--function betweenness` additionally needs the in-house igraph build in `vendor/` - see
[`vendor/README.md`](vendor/README.md) and the section below. That build is a superset
of stock igraph, so one environment can run every function.

**Import igraph before numpy.** On some clusters the in-house build was made with a newer
compiler than the system's and finds its newer `libstdc++` through its own library path.
If numpy loads the system copy first, `import igraph` fails with `GLIBCXX_3.4.32 not
found`. The package imports igraph first, so the command line and its worker processes
are safe; a script of your own that imports `numpy` before `graph_analysis` needs
`import igraph` first, or the newer `libstdc++` directory on `LD_LIBRARY_PATH`. On such a
cluster, load the Python module (`module load python/...`) before running, in batch scripts
and interactive sessions alike.

## Quick start

```bash
python3 src/main.py --function topology     --input-file examples/data/atrp_s1_sys1_0.data --output-file results.pkl
python3 src/main.py --function shortest-path --input-file examples/data/atrp_s1_sys1_0.data --ncpus 8 --output-file sp.pkl
python3 src/main.py --function betweenness  --input-file examples/data/atrp_s1_sys1_0.data --direction x --validate-stock
```

`--input-file` is a full LAMMPS data file (`write_data` output). Topology reads only its
Bonds section; shortest-path and betweenness also read the atom coordinates and the
orthorhombic box. `python3 src/main.py --help` lists every option.

As a library:

```python
import igraph  # first, see Installation
from graph_analysis import PolymerNetwork, NetworkAnalysisPipeline
from graph_analysis.shortest_path import PeriodicShortestPathAnalysis, ShortestPathConfig

network = PolymerNetwork.from_lammps("atrp_s1_sys1_0.data")
result = NetworkAnalysisPipeline(network).run()      # a ScanlanCaseResult: result.refined, result.raw

config = ShortestPathConfig(direction="x", max_period=4, ncpus=8)
sp = PeriodicShortestPathAnalysis.from_lammps("atrp_s1_sys1_0.data", config).run()
print(sp.min_length, len(sp.paths))

from graph_analysis.betweenness import SpatialBetweennessAnalysis  # needs the in-house build
bc = SpatialBetweennessAnalysis.from_lammps("atrp_s1_sys1_0.data").run()
```

See [`src/README.md`](src/README.md) for the full CLI reference and complete public API.

## Package layout

```
src/graph_analysis/
    network.py        PolymerNetwork - a networkx graph of the topology (all formats)
    lammps.py         LammpsDataReader, LammpsNetwork (atoms, coordinates, box), Direction
    crosslink.py      CrosslinkPathAnalyzer - strands and danglings by walking chains
    loops.py          LoopAnalyzer - minimum cycle basis loop detection
    scanlan_case.py   ScanlanCaseAnalyzer - crosslinker density refinement
    pipeline.py       NetworkAnalysisPipeline - runs all of the above in one pass
    shortest_path.py  PeriodicShortestPathAnalysis - shortest paths across the boundary
    betweenness.py    SpatialBetweennessAnalysis - native GEBC/OBC/DBC (needs vendor/)
    cli.py            command-line entry point (--function topology/shortest-path/betweenness)
examples/           data/ plus example_1, example_2, example_3 (see below)
vendor/
    igraph/           git submodule: in-house igraph C core (spatial-betweenness-1.0.0)
    python-igraph/    git submodule: in-house Python bindings (same branch)
```

Every graph view used by `--function topology` (the `networkx.Graph` used for
strand/dangling search and degree lookups, and the `igraph.Graph` used for loop
detection) is built from the same single, file-order-preserving parse in
`PolymerNetwork` - so node IDs (raw atom IDs from the topology file) are guaranteed
consistent across strands, danglings, and loops, and a loop's member atoms can be
cross-referenced against strand/dangling endpoints directly by ID. `lammps.py` is the
separate, array-based view that carries coordinates and the box for the spatial analyses;
the two are not interchangeable.

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
PyPI `igraph` doesn't have `edge_betweenness_spatial`; calling `--function betweenness` with it raises a clear
`RuntimeError` rather than silently falling back to something else.

### 3. Inspect new input before computing anything

```bash
vendor/.spatial-venv/bin/python src/main.py --function betweenness \
    --input-file examples/data/atrp_s1_sys1_0.data --inspect
```

Parses and validates the file, prints its topology/box as JSON, and does no
centrality calculation. Confirms atom/bond counts and the *actual* box lengths before
you commit to a multi-minute run - OBC and DBC depend on those lengths.

### 4. Validate on first use with new data, then drop `--validate-stock` for routine runs

```bash
# First run on a new network: also runs stock GEBC and fails loudly if native GEBC disagrees
vendor/.spatial-venv/bin/python src/main.py --function betweenness \
    --input-file examples/data/atrp_s1_sys1_0.data --direction x --validate-stock
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
with open("examples/example_2/atrp_s1_sys1_0.x.spatial.pkl", "rb") as f:
    result = pickle.load(f)  # needs graph_analysis importable, e.g. PYTHONPATH=src, with igraph imported first
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

## Examples

Each example folder has a `TUTORIAL.md` - a step-by-step walkthrough of one function on
one network - and the outputs it produces. All three use the same file,
`examples/data/atrp_s1_sys1_0.data` (103,149 atoms, 105,869 bonds, an ATRP network after
removing unreacted monomer at 96% conversion; the same network as the cluster's
`sys1/s1`).

| Folder | `--function` | Contents |
|---|---|---|
| [`examples/example_1/`](examples/example_1/TUTORIAL.md) | `topology` | strands, danglings, loops and the Scanlan-Case refinement |
| [`examples/example_2/`](examples/example_2/TUTORIAL.md) | `betweenness` | GEBC, OBC and DBC of every bond, with the check of the native GEBC against stock igraph |
| [`examples/example_3/`](examples/example_3/TUTORIAL.md) | `shortest-path` | the shortest paths across the periodic boundary, with the results of the cluster run they were checked against |

## Reproducibility notes

- **Node IDs are raw atom IDs everywhere** - the topology file's own numbering, used
  consistently as `networkx` node labels, `igraph` vertex indices, and every strand/
  dangling/loop/shortest-path entry. No remapping step is needed to cross-reference
  between them.
- **Loop identity depends on edge order, not vertex count.** A network's minimum cycle
  basis isn't always unique; when equally-short candidate cycles tie, `igraph`'s
  tie-breaking depends on the order edges were added, not on how many unused vertex IDs
  the graph was padded out to (`--nvertices`, if given, is otherwise safe to omit or
  over-specify). See `src/graph_analysis/loops.py` for the measurements behind this.
  Likewise the shortest-path graph is built in the reference code's edge order, since
  which of several equally short paths is returned depends on it.
- **Shortest paths reproduce the reference code.** For a given network the lengths and
  paths equal those of the Jackson lab code as run on the cluster, with the same
  `buffer_images`; the differences are that paths are reported in atom IDs and in images
  counted from image 0, and that atoms written unwrapped are handled (the reference
  silently drops their periodic bonds). See the module docstring of
  `src/graph_analysis/shortest_path.py`.
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

If you use `--function shortest-path`, please also cite the method it implements:

> Yu, Z.; Jackson, N. E. Shortest Paths Govern Bond Rupture in Thermoset Networks.
> *Macromolecules* **2025**.
> DOI: [10.1021/acs.macromol.4c02438](https://doi.org/10.1021/acs.macromol.4c02438)

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

@article{yu2025shortest,
  title   = {Shortest Paths Govern Bond Rupture in Thermoset Networks},
  author  = {Yu, Zheng and Jackson, Nicholas E.},
  journal = {Macromolecules},
  year    = {2025},
  doi     = {10.1021/acs.macromol.4c02438}
}
```
