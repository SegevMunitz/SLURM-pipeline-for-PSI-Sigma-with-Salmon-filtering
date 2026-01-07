#!/usr/bin/env python3
"""PSI-Sigma runner (Python wrapper around dummyai.pl).

Run on an HPC node (via SLURM). This script:
  - validates inputs (GTF, group lists, PSI-Sigma script)
  - configures a user-local Perl installation (PERL5LIB/PATH)
  - creates a per-job TMPDIR
  - removes existing output directories (PSI-Sigma fails if they already exist)
  - runs one or more PSI-Sigma comparisons

Edit paths.py to match your environment.
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
    if not p.is_file():
        raise SystemExit(f"ERROR: missing file: {p}")


def build_env(perlbase: Path) -> dict[str, str]:
    """Return a copy of os.environ with PERL5LIB/PATH updated."""
    env = dict(os.environ)

    perl_lib = perlbase / "lib" / "perl5"
    perl_bin = perlbase / "bin"

    env["PERL5LIB"] = f"{perl_lib}:{env.get('PERL5LIB', '')}".rstrip(":")
    env["PATH"] = f"{perl_bin}:{env.get('PATH', '')}"

    return env


def remove_if_exists(p: Path) -> None:
    if p.exists():
        print(f"Removing existing output dir: {p}")
        shutil.rmtree(p)


def main() -> None:
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
