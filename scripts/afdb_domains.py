#!/usr/bin/env python3
"""
Executable shell for the afdb-domains CLI.
==========================================

Adds the repository root to ``sys.path`` so the flat ``src`` package can be
imported when the script is invoked by path. All argument parsing, validation,
and error handling lives in :mod:`src.cli`.

Usage
-----
    python scripts/afdb_domains.py --help
    python scripts/afdb_domains.py pae-to-domains --pae-dir <dir> --output-dir <dir>

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-11-03
Version: 2.0.0
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.cli import main

# =============================================================================
# 2. ENTRY POINT
# =============================================================================
if __name__ == "__main__":
    sys.exit(main())
