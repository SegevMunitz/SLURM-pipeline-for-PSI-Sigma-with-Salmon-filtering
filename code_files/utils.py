from __future__ import annotations
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def which_or_die(exe: str) -> str:
    p = shutil.which(exe)
    if not p:
        raise SystemExit(f"ERROR: '{exe}' not found in PATH.")
    return p

def run(cmd: List[str], cwd: Optional[Path] = None, env: Optional[dict[str, str]] = None) -> None:
    print("\n$ " + " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, env=env, check=True)

def write_text(p: Path, s: str) -> None:
    ensure_dir(p.parent)
    p.write_text(s, encoding="utf-8")
