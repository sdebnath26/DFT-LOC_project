#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from pathlib import Path

# NOTE:
# This module was originally exported from a Jupyter notebook.
# In CI (GitHub Actions), IPython/Jupyter helpers like `get_ipython()` are not available.
# Keep this file import-safe by guarding such calls.
try:
    from IPython import get_ipython  # type: ignore


  _ip = get_ipython()
    if _ip is not None:
        _ip.run_line_magic("load_ext", "autoreload")
        _ip.run_line_magic("autoreload", "2")
except Exception:
    # Not running inside IPython/Jupyter.
    pass

file_path = Path(r"Propane-new.xyz")

def open_xyz(path: str | Path, encoding: str = "utf-8") -> tuple[list[str], list[str]]:
    """Read an .xyz file and return (rows, data)."""
