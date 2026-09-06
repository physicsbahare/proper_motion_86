"""Compatibility shim for legacy first26 imports.

The validated first26 recentering helper historically imports
``forced_astrometry`` as a top-level module.  Re-export the canonical
implementation so package-style execution (``python -m defensible.run_catalog``)
works in GitHub Actions without changing the scientific code path.
"""
from first26.forced_astrometry import *  # noqa: F401,F403
