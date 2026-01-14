#!/usr/bin/env python3
"""
step5_run_psi_sigma.py

PSI-Sigma runner for alternative splicing / junction usage comparisons (Perl wrapper).

This script is a thin, SLURM-friendly Python wrapper around the PSI-Sigma Perl
pipeline (dummyai.pl). It is intended to run on an HPC compute node after STAR
alignment has completed and after group files have been generated. Group files
are plain text lists containing one absolute BAM path per line.

The script performs four key roles:
  1) Input validation: verifies required reference files (GTF) and all group lists
     exist, and confirms required executables (perl) are available.
  2) Environment configuration: injects a user-local Perl installation into PATH
     and PERL5LIB (via paths.PERLBASE) so PSI-Sigma can run without system-wide
     Perl modules.
  3) Execution hygiene: creates a per-job TMPDIR to avoid collisions, ensures the
     pipeline output/log directories exist, and removes PSI-Sigma output folders
     if they already exist (PSI-Sigma commonly fails when outputs pre-exist).
  4) Batch execution: runs one or more comparisons defined in paths.COMPARISONS,
     each producing a separate output directory under paths.OUTDIR/<comparison>.

Configuration is controlled entirely via paths.py:
  - ANNOTATION_GTF: reference annotation used by PSI-Sigma
  - PSI_SIGMA: path to dummyai.pl (PSI-Sigma entry point)
  - PERLBASE: root of user-local Perl installation (bin/, lib/perl5/)
  - OUTDIR: pipeline output root
  - COMPARISONS: list of comparison specs, each with:
        name, groupa (Path), groupb (Path)

Notes
-----
- Threads are taken from SLURM_CPUS_PER_TASK when available, otherwise default to 4.
- TMPDIR is set per job to avoid shared /tmp collisions on compute nodes.
- This script intentionally removes existing comparison output directories to
  prevent PSI-Sigma from failing due to pre-existing outputs.

Typical usage
-------------
Run as a single SLURM job after alignment + group file generation, e.g.:
  sbatch sbatch_05_run_psi_sigma.slurm
"""

from __future__ import annotations

import os
import shutil
import socket
from pathlib import Path

import paths
from utils import ensure_dir, run, which_or_die


# =========================
# PSI-Sigma parameters (tune as needed)
# =========================
PSI_TYPE = "1"
NREAD = "10"
IRRANGE = "0"
FMODE = "3"


def ensure_file(p: Path) -> None:
    """
    Validate that a required input file exists on disk.

    Used for reference inputs and group lists to fail fast with a clear error
    message before starting long-running PSI-Sigma computations.
    """
    if not p.is_file():
        raise SystemExit(f"ERROR: missing file: {p}")


def build_env(perlbase: Path) -> dict[str, str]:
    """
    Construct an execution environment for PSI-Sigma with a user-local Perl.

    Returns a copy of os.environ with PERL5LIB updated to include perlbase/lib/perl5
    and PATH updated to include perlbase/bin, enabling PSI-Sigma dependencies to
    resolve without relying on system Perl modules.
    """
    env = dict(os.environ)

    perl_lib = perlbase / "lib" / "perl5"
    perl_bin = perlbase / "bin"

    env["PERL5LIB"] = f"{perl_lib}:{env.get('PERL5LIB', '')}".rstrip(":")
    env["PATH"] = f"{perl_bin}:{env.get('PATH', '')}"

    return env


def remove_if_exists(p: Path) -> None:
    """
    Remove an existing directory tree if present.

    PSI-Sigma often fails if the target output directory already exists, so this
    helper enforces a clean output location for each comparison run.
    """
    if p.exists():
        print(f"Removing existing output dir: {p}")
        shutil.rmtree(p)


def main() -> None:
    """
    Validate configuration, set up environment/TMPDIR, and run PSI-Sigma comparisons.

    Uses paths.py to locate GTF, PSI-Sigma script, Perl base, and group files.
    For each comparison in paths.COMPARISONS, removes any pre-existing output
    directory and launches dummyai.pl with consistent parameters and threading.
    """
    # Threads: prefer SLURM, fallback to 4
    n_threads = os.environ.get("SLURM_CPUS_PER_TASK", "4")
    job_id = os.environ.get("SLURM_JOB_ID", "local")

    print(f"=== job_id={job_id} host={socket.gethostname()} ===")
    print(f"Threads: {n_threads}")

    # Validate required files
    ensure_file(paths.ANNOTATION_GTF)
    ensure_file(paths.PSI_SIGMA)

    for comp in paths.COMPARISONS:
        ensure_file(comp.groupa)
        ensure_file(comp.groupb)

    which_or_die("perl")

    # Environment setup (Perl)
    env = build_env(paths.PERLBASE)

    # TMPDIR: avoid shared /tmp collisions
    tmpdir = Path(env.get("TMPDIR", str(paths.OUTDIR / "tmp" / f"run_{job_id}")))
    ensure_dir(tmpdir)
    env["TMPDIR"] = str(tmpdir)

    # Output dirs
    ensure_dir(paths.OUTDIR)
    ensure_dir(paths.OUTDIR / "logs")

    # PSI-Sigma fails if output directories already exist
    for comp in paths.COMPARISONS:
        remove_if_exists(paths.OUTDIR / comp.name)

    print("PSI-Sigma script:", paths.PSI_SIGMA)
    print("GTF:", paths.ANNOTATION_GTF)
    print("OUT_DIR:", paths.OUTDIR)
    print("TMPDIR:", tmpdir)

    # Run comparisons
    for comp in paths.COMPARISONS:
        print(f"\n=== Running PSI-Sigma: {comp.name} ===")

        out = paths.OUTDIR / comp.name
        cmd = [
            "perl", str(paths.PSI_SIGMA),
            "--gtf", str(paths.ANNOTATION_GTF),
            "--name", comp.name,
            "--type", PSI_TYPE,
            "--nread", NREAD,
            "--irrange", IRRANGE,
            "--fmode", FMODE,
            "--groupa", str(comp.groupa),
            "--groupb", str(comp.groupb),
            "--output", str(out),
            "--threads", str(n_threads),
        ]
        run(cmd, cwd=None, env=env)

    print("\n=== DONE ===")


if __name__ == "__main__":
    main()
