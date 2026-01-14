#!/usr/bin/env python3
"""
utils.py

Shared utility helpers for the RNA-seq SLURM pipeline.

This module provides small, reusable helper functions used across multiple
pipeline steps, including directory creation, executable validation, command
execution, and safe text-file writing. The goal is to centralize common
boilerplate and enforce consistent behavior (error handling, logging) across
all pipeline scripts.

These utilities are intentionally lightweight and dependency-free so they can
be imported safely in both interactive runs and SLURM batch jobs.
"""

from __future__ import annotations
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

def ensure_dir(p: Path) -> None:
    """
    Ensure that a directory exists.

    Creates the directory and any missing parent directories if needed.
    If the directory already exists, no action is taken.
    """
    p.mkdir(parents=True, exist_ok=True)

def which_or_die(exe: str) -> str:
    """
    Locate an executable in PATH or terminate with a clear error.

    Uses shutil.which to resolve the executable. If not found, exits
    immediately to prevent silent failures later in the pipeline.
    """
    p = shutil.which(exe)
    if not p:
        raise SystemExit(f"ERROR: '{exe}' not found in PATH.")
    return p

def run(cmd: List[str], cwd: Optional[Path] = None, env: Optional[dict[str, str]] = None) -> None:
    """
    Execute an external command with logging and strict error checking.

    Prints the command (shell-style) for transparency, then runs it using
    subprocess.run with check=True so failures propagate immediately.
    """
    print("\n$ " + " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, env=env, check=True)

def write_text(p: Path, s: str) -> None:
    """
    Write text content to a file, creating parent directories if needed.

    Ensures the parent directory exists before writing and encodes output
    as UTF-8 to avoid platform-dependent defaults.
    """
    ensure_dir(p.parent)
    p.write_text(s, encoding="utf-8")
