"""Import order, which decides whether igraph can find its C++ runtime on some clusters."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SRC = str(Path(__file__).resolve().parents[1] / "src")


def test_igraph_is_imported_before_numpy_and_networkx():
    code = (
        "import sys; sys.path.insert(0, %r)\n"
        "import graph_analysis\n"
        "order = list(sys.modules)\n"
        "assert order.index('igraph') < order.index('numpy'), 'numpy was imported before igraph'\n"
        "assert order.index('igraph') < order.index('networkx'), 'networkx was imported before igraph'\n"
    ) % SRC
    subprocess.run([sys.executable, "-I", "-c", code], check=True)
