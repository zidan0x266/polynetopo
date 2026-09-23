#!/usr/bin/env python3
"""Convenience entry point: `python3 main.py <command> ...` without installing the package.

See `graph_analysis.cli` for the actual implementation and `python3 main.py --help`
for the available subcommands.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from graph_analysis.cli import main

if __name__ == "__main__":
    main()
