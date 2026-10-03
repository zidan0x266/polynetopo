"""Native spatially-weighted edge betweenness (GEBC/OBC/DBC) for LAMMPS polymer networks.

Computes geodesic edge betweenness centrality (GEBC, the standard topological
metric) alongside two spatially-weighted variants in a single pass: OBC
(orientation-weighted, by alignment with a loading axis) and DBC
(distance-weighted, by source-target separation) - per Bai, Ding, Qi, Peng,
Huang, Yan & Man, "Spatially Weighted Betweenness Centrality for Predicting
Fracture in End-Linked Polymer Networks" (Macromolecules, 2026).

Requires the in-house igraph build with `Graph.edge_betweenness_spatial` (see
`vendor/README.md` at the repo root for how to build it) - the stock PyPI
igraph package does not have this method, and importing this subpackage
without it raises a clear `RuntimeError` at analysis time, not at import time.

This subpackage's `PolymerNetwork` is a distinct, array-based representation
(atom IDs, coordinates, box bounds) from `graph_analysis.network.PolymerNetwork`
(a `networkx`-backed graph for the topology/shortest-path analyses) - the two
aren't interchangeable, so it's deliberately not re-exported from the top-level
`graph_analysis` package. Import it from here explicitly:
`from graph_analysis.betweenness import SpatialBetweennessAnalysis`.
"""

from .analysis import CentralityValidationError, SpatialBetweennessAnalysis
from .lammps import LammpsDataReader, LammpsFormatError
from .models import AnalysisConfig, Direction, PolymerNetwork
from .results import CentralityResult, ValidationReport

__all__ = [
    "AnalysisConfig",
    "CentralityResult",
    "CentralityValidationError",
    "Direction",
    "LammpsDataReader",
    "LammpsFormatError",
    "PolymerNetwork",
    "SpatialBetweennessAnalysis",
    "ValidationReport",
]
