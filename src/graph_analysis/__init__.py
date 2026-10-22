"""Topology analysis pipeline for crosslinked polymer networks.

Public API:
    PolymerNetwork        - load a network topology (RMC or LAMMPS format)
    CrosslinkAction        - CROSSLINKERS or DANGLINGS search presets
    CrosslinkPathAnalyzer  - strands and dangling chains between crosslinkers
    LoopAnalyzer, LoopBasis - minimum cycle basis / loop detection
    NetworkStats           - strands/danglings/loops for one network snapshot
    RefinedStats            - NetworkStats plus the crosslinker density it produced
    ScanlanCaseAnalyzer, ScanlanCaseResult - crosslinker density refinement
    NetworkAnalysisPipeline - runs all of the above in one pass
"""

# igraph must be imported before numpy/networkx. A build made with a newer compiler than the
# system's needs the newer libstdc++ and finds it through its own library path; if numpy has
# already loaded the system copy, importing igraph fails with "GLIBCXX_3.4.32 not found" (seen on
# a cluster node). Worker processes import this package first, so they get the same order.
import igraph  # noqa: F401, I001

from .crosslink import CrosslinkAction, CrosslinkPathAnalyzer
from .loops import LoopAnalyzer, LoopBasis
from .network import PolymerNetwork
from .pipeline import NetworkAnalysisPipeline
from .scanlan_case import (
    NetworkStats,
    RefinedStats,
    ScanlanCaseAnalyzer,
    ScanlanCaseResult,
    merge_strands,
)

__all__ = [
    "PolymerNetwork",
    "CrosslinkAction",
    "CrosslinkPathAnalyzer",
    "LoopAnalyzer",
    "LoopBasis",
    "NetworkStats",
    "RefinedStats",
    "ScanlanCaseAnalyzer",
    "ScanlanCaseResult",
    "merge_strands",
    "NetworkAnalysisPipeline",
]
