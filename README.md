<p align="center">
  <img src="logo/polynetopo_logo.png?raw=true" alt="polynetopo" width="640">
</p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.24325718"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.24325718.svg" alt="DOI 10.5281/zenodo.24325718"></a>
</p>

# polynetopo

Topology, shortest-path and spatial betweenness analysis for crosslinked polymer
networks:

<table>
  <thead>
    <tr>
      <th align="left"><code>--function</code></th>
      <th align="left">What it computes</th>
      <th align="left">Needs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td nowrap><a href="#--function-topology"><b><samp>topology</samp></b></a></td>
      <td>Elastically effective strands, dangling chains, primary/higher-order loops, and the Scanlan-Case crosslinker density refinement</td>
      <td nowrap><img alt="bonds" src="https://img.shields.io/badge/bonds-2f81f7"></td>
    </tr>
    <tr>
      <td nowrap><a href="#--function-shortest-path"><b><samp>shortest-path</samp></b></a></td>
      <td>The shortest paths that wrap through the periodic box, which govern where and when the first bond breaks</td>
      <td nowrap><img alt="bonds" src="https://img.shields.io/badge/bonds-2f81f7"><br><img alt="coordinates" src="https://img.shields.io/badge/coordinates-2da44e"><br><img alt="box" src="https://img.shields.io/badge/box-8250df"></td>
    </tr>
    <tr>
      <td nowrap><a href="#--function-betweenness"><b><samp>betweenness</samp></b></a></td>
      <td>Native GEBC / OBC / DBC edge betweenness</td>
      <td nowrap><img alt="bonds" src="https://img.shields.io/badge/bonds-2f81f7"><br><img alt="coordinates" src="https://img.shields.io/badge/coordinates-2da44e"><br><img alt="box" src="https://img.shields.io/badge/box-8250df"></td>
    </tr>
  </tbody>
</table>

This code produced the analysis in the three papers listed first under
[Citing this work](#citing-this-work). It is built on `networkx` and `igraph`; the betweenness analysis needs the in-house
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
- **Scanlan-Case-refined crosslinker density** - keeps the load-bearing part of the
  network, its largest biconnected component, without everything that hangs on it through
  a single strand or a single crosslinker (dangling chains, primary loops and dangling
  clusters). A crosslinker that carries a primary loop or dangling chain(s) thus becomes an
  f=3 junction or a linear unit, the strands meeting at a linear unit are spliced into one,
  the strands that lead only to a dead end join the dangling chains, and the effective
  vs. total crosslinker counts are recomputed.

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

A shortest path (SP) that wraps through the periodic box can be longer than one image,
so for each period p = 1, 2, ..., M (`--max-period`, default 8):

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
   of a path's atoms. This covers Case 1 (the same path found from
   different start atoms) and, together with the validity rule, Case 2 (a path at larger
   p that repeats one at smaller p). `--keep-redundant` skips it.

The search stops at the first period with no valid path. The smallest D^g is the first
(shortest) SP; it predicts where the first bond breaks and, through
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
must have no repeated bonds (an atom pair bonded twice); this raises an error naming the
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

## Installation

Requires Python >= 3.10 with `networkx` (>= 3.0), `numpy`, `igraph` (the `python-igraph`
package) and `tqdm`. Install with `pip install -e .`, which provides the `polynetopo`
command (the Python package is imported as `graph_analysis`), or skip the install and run
`python3 src/main.py` from the repository root, which works without it.
`--function betweenness` additionally needs the in-house igraph build in `vendor/` - see
[`vendor/README.md`](vendor/README.md) and the section below. That build is a superset
of stock igraph, so one environment can run every function.

## Tests

```bash
pip install -e ".[test]"
python3 -m pytest
```

The tests take about a second. Those of `--function betweenness` need the in-house igraph
build and are skipped without it; to run them, use that build's Python
(`vendor/.spatial-venv/bin/python -m pip install pytest`, then
`vendor/.spatial-venv/bin/python -m pytest`).

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
tests/              pytest suite (see Tests)
vendor/
    igraph/           git submodule: in-house igraph C core (spatial-betweenness-1.0.0)
    python-igraph/    git submodule: in-house Python bindings (same branch)
```

## Using `--function betweenness` correctly

This section covers the Python interpreter to use, the validation on new data, and
the raw DBC convention.

### 1. Build once per machine

```bash
git submodule update --init vendor/igraph vendor/python-igraph
echo "1.0.0" > vendor/igraph/IGRAPH_VERSION   # the version string for the build; see vendor/README.md
cd vendor && bash python-igraph/scripts/build_spatial.sh && cd ..
vendor/.spatial-venv/bin/python -m pip install --no-build-isolation --no-deps networkx tqdm
vendor/.spatial-venv/bin/python -m pip install -e .
```

Takes ~15-20 minutes (full C+Python build and test suite); do this once per machine -
binaries from one machine/OS don't work on another. Full details are in
[`vendor/README.md`](vendor/README.md).

### 2. Always use that build's Python, never your regular one

```bash
vendor/.spatial-venv/bin/python -c 'import igraph as ig; print(hasattr(ig.Graph, "edge_betweenness_spatial"))'
# must print True
```

Every command below uses `vendor/.spatial-venv/bin/python`, not plain `python3`; stock
PyPI `igraph` doesn't have `edge_betweenness_spatial`.

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

Output is a pickled `CentralityResult`:

```python
import pickle
with open("examples/example_2/atrp_s1_sys1_0.x.spatial.pkl", "rb") as f:
    result = pickle.load(f)  # needs graph_analysis importable, e.g. PYTHONPATH=src
raw_gebc, raw_obc, raw_dbc = result.raw_gebc, result.raw_obc, result.raw_dbc
obc, dbc = result.obc, result.dbc  # pair-normalized properties; dbc also length-divided IF --dbc-length was passed
```

`raw_dbc` is **always** the distance-only sum (no `Lx`/bond-length denominator) -
apply any reference length yourself, as many times as you like, without recomputing:
`result.dbc_for_length(L)` or `result.dbc_for_length(L, normalized=True)`. Passing
`--dbc-length L` at the CLI only changes the *postprocessed* `dbc` property; `raw_dbc`
is unaffected either way. See `vendor/python-igraph/analysis/README.md` for the full
convention reference (normalization factors, what each field indexes).

OBC and DBC weight each source-target pair by its minimum-image displacement, whereas the
shortest path between the two atoms on the periodic graph may end on another image of the
target; this happens more often the longer the path is compared with the box.

## Examples

Each example folder has a `TUTORIAL.md` - a step-by-step walkthrough of one function on
one network - and the outputs it produces. All three use the same file,
`examples/data/atrp_s1_sys1_0.data` (103,149 atoms, 105,869 bonds, an ATRP network after
removing unreacted monomer at 96% conversion).

| Folder | `--function` | Contents |
|---|---|---|
| [`examples/example_1/`](examples/example_1/TUTORIAL.md) | `topology` | strands, danglings, loops and the Scanlan-Case refinement |
| [`examples/example_2/`](examples/example_2/TUTORIAL.md) | `betweenness` | GEBC, OBC and DBC of every bond, with the check of the native GEBC against stock igraph |
| [`examples/example_3/`](examples/example_3/TUTORIAL.md) | `shortest-path` | the shortest paths across the periodic boundary of the network |

## Citing this work

If you use this code, please cite:

> Zhang, Z.; Krajniak, J.; Dookhith, A. Z.; Tian, Y.; Sachar, H. S.; Marioni, N.;
> Duncan, T. J.; Liu, J.; Sanoja, G. E.; Ganesan, V. Topology and Mechanical Properties
> of Polymer Networks Formed under Free Radical and Atom Transfer Radical
> Polymerizations. *Macromolecules* **2025**, *58*, 3168-3187.
> DOI: [10.1021/acs.macromol.4c03138](https://doi.org/10.1021/acs.macromol.4c03138)

> Zhang, Z.; Krajniak, J.; Dookhith, A. Z.; Zhang, H.; Tian, Y.; Sachar, H. S.; Marioni,
> N.; Duncan, T. J.; Liu, J.; Sanoja, G. E.; Ganesan, V. Fracture Behavior of
> Monodisperse and Polydisperse Polymer Networks: Effects of Topology and
> Strand-Length Dispersity. *Macromolecules* **2026**, *59*, 99-119.
> DOI: [10.1021/acs.macromol.5c01715](https://doi.org/10.1021/acs.macromol.5c01715)

> Zhang, Z.; Krajniak, J.; Tian, Y.; Wang, S.; Mukherjee, K.; Liu, J.; Chen, W.; Ganesan, V.
> Are Short Strands the Weak Links of Polymer Networks? Taut Paths and the Fracture Energy
> of Polydisperse Networks. In preparation.

The software itself (GitHub's "Cite this repository" gives this entry in APA and BibTeX):

> Zhang, Z.; Krajniak, J.; Wang, S.; Ganesan, V. polynetopo: Topology, Shortest-Path and
> Spatial Betweenness Analysis for Crosslinked Polymer Networks, version 0.5.0. Zenodo,
> **2026**. DOI: [10.5281/zenodo.24325718](https://doi.org/10.5281/zenodo.24325718)

Optionally, depending on the functions you use:

If you use the shortest paths (`--function shortest-path`), please also cite the method it
implements:

> Yu, Z.; Jackson, N. E. Shortest Paths Govern Bond Rupture in Thermoset Networks.
> *Macromolecules* **2025**, *58*, 1728-1736.
> DOI: [10.1021/acs.macromol.4c02438](https://doi.org/10.1021/acs.macromol.4c02438)

If you use the weighted betweenness (`--function betweenness`), please also cite:

> Bai, X.; Ding, Z.; Qi, S.; Peng, L.; Huang, X.; Yan, D.; Man, X. Spatially Weighted
> Betweenness Centrality for Predicting Fracture in End-Linked Polymer Networks.
> *Macromolecules* **2026**, *59*, 10294-10303.
> DOI: [10.1021/acs.macromol.6c00974](https://doi.org/10.1021/acs.macromol.6c00974)

## License

polynetopo is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version; see [`LICENSE`](LICENSE).
The in-house igraph build in `vendor/` keeps the license of igraph, the GNU General
Public License, version 2 or later.
