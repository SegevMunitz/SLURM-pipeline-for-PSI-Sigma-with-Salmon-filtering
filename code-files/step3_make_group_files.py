#!/usr/bin/env python3
"""
step3_make_group_files.py

Create PSI-Sigma group files (lists of BAM paths) after STAR alignment.

What this script does
---------------------
Given:
  1) an output root directory (OUTDIR) that contains STAR outputs in:
       <OUTDIR>/bam/<sample>/<sample>.Aligned.sortedByCoord.out.bam
  2) a samples file:
       <OUTDIR>/samples.txt
     containing one sample ID per line (the same IDs used in paths.SAMPLES / step1 prep)
  3) a mapping TSV file:
       groups.tsv
     with exactly 2 tab-separated columns per line:
       <sample_id>    <group_name>

This script writes, for each group_name, a file:
  <OUTDIR>/groups/<group_name>.bams.txt

Each output file contains absolute BAM file paths, one per line. These files can be
passed directly to PSI-Sigma as group inputs (e.g., --groupa /path/to/Healthy.bams.txt).

Important usage note
--------------------
Run this as a SINGLE job after the STAR array has completed successfully.
Do NOT run it inside the STAR per-sample array task, because you can get partial
group files while alignment is still running.

Exit conditions / validation
----------------------------
- If a mapping line does not have exactly 2 TSV columns -> exits with error.
- If a mapped sample has no BAM file at the expected location -> exits with error.
- If a mapped sample is not present in samples.txt -> prints WARNING to stderr
  (but still includes it in the group; BAM existence is still enforced).
"""

from __future__ import annotations

from pathlib import Path
import argparse
import sys
from paths import (OUTDIR)

def bam_path(out_dir: Path, sample: str) -> Path:
    """
    Construct the expected STAR BAM path for a given sample.

    Returns the canonical location under the pipeline output directory:
      <out_dir>/bam/<sample>/<sample>.Aligned.sortedByCoord.out.bam
    """
    return out_dir / "bam" / sample / f"{sample}.Aligned.sortedByCoord.out.bam"


def main() -> None:
    """
    Create group-specific BAM list files from sample-to-group mappings.

    Reads sample IDs and group assignments, validates BAM existence, and
    writes one BAM list file per group under <OUTDIR>/groups for downstream
    PSI-Sigma analysis.
    """

    out_dir = OUTDIR / "star_alignments"
    samples_txt = OUTDIR / "samples.txt" 
    groups_tsv = out_dir / "groups.tsv"
    groups_dir = out_dir / "groups"
    groups_dir.mkdir(parents=True, exist_ok=True)

    # read samples
    samples = [ln.strip() for ln in samples_txt.read_text().splitlines() if ln.strip()]
    sample_set = set(samples)

    # read mapping
    mapping_lines = [ln.strip() for ln in groups_tsv.read_text().splitlines() if ln.strip()]
    group_to_samples: dict[str, list[str]] = {}
    for ln in mapping_lines:
        parts = ln.split("\t")
        if len(parts) != 2:
            raise SystemExit(f"Bad mapping line (need 2 TSV columns): {ln}")
        s, g = parts[0].strip(), parts[1].strip()
        if s not in sample_set:
            print(f"WARNING: mapping sample not in samples.txt: {s}", file=sys.stderr)
        group_to_samples.setdefault(g, []).append(s)

    # write group files
    written = []
    for group, ss in group_to_samples.items():
        p = groups_dir / f"{group}.bams.txt"
        bams = []
        for s in ss:
            b = bam_path(out_dir, s)
            if not b.is_file():
                raise SystemExit(f"Missing BAM for group={group} sample={s}: {b}")
            bams.append(str(b.resolve()))
        p.write_text("\n".join(bams) + "\n")
        written.append(p)

    print("Wrote group files:")
    for p in written:
        print("  ", p)


if __name__ == "__main__":
    main()
