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
