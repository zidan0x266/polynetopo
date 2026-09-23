# PolyNetopo

Topology analysis for crosslinked polymer networks: elastically effective strands,
dangling chains, primary/higher-order loops, and the Scanlan-Case crosslinker density
refinement. Built on `networkx` (shortest-path search) and `igraph` (minimum cycle
basis), with a multiprocessing-accelerated search over large networks.

This code produced the topological analysis in the two papers cited at the bottom of
this README.

## What it does

Given a polymer network topology (a bond/edge list - from a LAMMPS bond table, an
RMC-style edge list, or a reactive Monte Carlo output converted to one of those), this
package identifies:

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

## Installation

```bash
pip install -e .
```

Requires Python >= 3.10, `networkx`, `numpy`, `igraph` (the `python-igraph` package),
and `tqdm`.

## Quick start

Run the whole pipeline - strands, danglings, loops, and their Scanlan-Case refinement -
from a single topology load, with no intermediate files:

```bash
python3 src/main.py analyze --input-file examples/atrp/edges_0.bds \
    --output-file results.pkl --raw --ncpus 4 --max-chunk-size 10000
```

This writes `results.pkl` (a `graph_analysis.RefinedStats`: the refined
strands/danglings/loops plus the crosslinker density they produced) and, since `--raw`
was passed, `results_raw.pkl` (a `graph_analysis.NetworkStats`: the pre-refinement data,
for sanity-checking a configuration without redoing the analysis).

As a library:

```python
from graph_analysis import PolymerNetwork, NetworkAnalysisPipeline, PipelineConfig

network = PolymerNetwork.from_lammps("edges_0.bds")
pipeline = NetworkAnalysisPipeline(network, PipelineConfig(sys_name="run1", ncpus=4))
result = pipeline.run()  # a ScanlanCaseResult: result.refined, result.raw
```

See [`src/README.md`](src/README.md) for the full CLI reference (including the
step-by-step subcommands - `distances`, `loops`, `scanlan-case`, `compare` - useful for
inspecting an intermediate stage or resuming after a crash) and the complete public API.

## Package layout

```
src/graph_analysis/
    network.py       PolymerNetwork - loads a topology (RMC or LAMMPS format)
    crosslink.py      CrosslinkPathAnalyzer - shortest simple paths between crosslinkers
    loops.py          LoopAnalyzer - minimum cycle basis loop detection
    scanlan_case.py   ScanlanCaseAnalyzer - crosslinker density refinement
    pipeline.py       NetworkAnalysisPipeline - runs all of the above in one pass
    cli.py            command-line entry point
```

Every graph view (the `networkx.Graph` used for strand/dangling search and degree
lookups, and the `igraph.Graph` used for loop detection) is built from the same
single, file-order-preserving parse in `PolymerNetwork` - so node IDs (raw atom IDs
from the topology file) are guaranteed consistent across strands, danglings, and loops,
and a loop's member atoms can be cross-referenced against strand/dangling endpoints
directly by ID.

## Example dataset

`examples/atrp/` contains one worked example: an ATRP-synthesized network
(`edges_0.bds`, a LAMMPS bond table, 103,149 atoms after removing unreacted monomer at
96% conversion) and its analysis outputs (`strands.pkl`, `danglings.pkl`, `igloops.pkl`,
`Scanlan_Case.pkl`), reproducible with the commands in `src/README.md`.

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
```
