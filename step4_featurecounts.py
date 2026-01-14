#!/usr/bin/env python3
"""
step4_featurecounts.py

Run featureCounts on STAR-aligned BAMs and produce a PSI-Sigma-friendly count matrix.

This step aggregates per-sample alignments into gene-level read counts using
Subread's featureCounts. It is intended to run after alignment has completed
and the output directory contains:
    <OUT_ROOT>/bam/<sample>/<sample>.Aligned.sortedByCoord.out.bam

The script can either:
  - auto-discover BAMs under <OUT_ROOT>/bam/*/*.Aligned.sortedByCoord.out.bam, or
  - accept one or more BAMs explicitly via repeated --bam arguments.

Outputs are written under:
    <OUT_ROOT>/counts/
including:
  - featureCounts.gene_counts.txt   (raw featureCounts output)
  - gene_counts.matrix.tsv          (cleaned matrix)

The "clean" matrix is produced by removing featureCounts metadata columns
(Chr, Start, End, Strand, Length) and renaming the Geneid column to gene_id,
while preserving the sample count columns for downstream differential analysis.

Requirements
------------
- featureCounts must be available on PATH (module load / conda environment).
- The provided --annotation-gtf should match the genome used for alignment.
"""

from __future__ import annotations
import argparse
import csv
from pathlib import Path
from typing import List
from paths import (OUTDIR, ANNOTATION_GTF)
from utils import ensure_dir, run, which_or_die

def featurecounts_gene_counts(
    bams: List[Path],
    annotation_gtf: Path,
    out_root: Path,
    threads: int,
    stranded: int,
    paired_end: bool,
) -> Path:
    """
    Run featureCounts on the provided BAMs and write outputs under out_root/counts.

    Produces the raw featureCounts output file and then generates a cleaned
    tab-delimited matrix (gene_counts.matrix.tsv) suitable for PSI-Sigma.
    Returns the path to the cleaned matrix file.
    """
    which_or_die("featureCounts")
    out_dir = out_root / "counts"
    ensure_dir(out_dir)

    fc_out = out_dir / "featureCounts.gene_counts.txt"
    cmd = [
        "featureCounts",
        "-T", str(threads),
        "-a", str(annotation_gtf),
        "-o", str(fc_out),
        "-g", "gene_id",
        "-t", "exon",
        "-s", str(stranded),
    ]
    if paired_end:
        cmd += ["-p", "-B", "-C"]
    cmd += [str(b) for b in bams]
    run(cmd)

    # Create a clean matrix TSV for Psi-Sigma
    matrix = out_dir / "gene_counts.matrix.tsv"
    make_clean_counts_matrix(fc_out, matrix)
    return matrix

def make_clean_counts_matrix(featurecounts_out: Path, clean_out: Path) -> None:
    """
    Convert raw featureCounts output into a compact gene-by-sample count matrix.

    Skips comment lines, identifies the header row, drops featureCounts metadata
    columns (Chr/Start/End/Strand/Length), and renames Geneid -> gene_id.
    Writes the resulting TSV matrix to clean_out.
    """
    with open(featurecounts_out, "r", encoding="utf-8") as f_in, \
         open(clean_out, "w", encoding="utf-8", newline="") as f_out:
        
        reader = csv.reader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")

        # Skip comment lines (starting with #) until we find the header
        headers = []
        for row in reader:
            if row and not row[0].startswith("#"):
                headers = row
                break
        
        if not headers:
            raise SystemExit(f"ERROR: Could not find header row in {featurecounts_out}")

        # Identify columns to keep (Geneid + samples) and rename Geneid -> gene_id
        exclude = {"Chr", "Start", "End", "Strand", "Length"}
        keep_idxs = [i for i, h in enumerate(headers) if h == "Geneid" or h not in exclude]
        new_headers = [("gene_id" if headers[i] == "Geneid" else headers[i]) for i in keep_idxs]
        
        writer.writerow(new_headers)
        for row in reader:
            # Be defensive: skip malformed lines
            if len(row) < max(keep_idxs) + 1:
                continue
            writer.writerow([row[i] for i in keep_idxs])


def find_star_bams(out_root: Path) -> List[Path]:
    """
    Discover STAR-sorted BAMs under the pipeline output directory.

    Searches for files matching:
      <out_root>/bam/*/*.Aligned.sortedByCoord.out.bam
    Returns a sorted list of Paths for reproducible ordering.
    """
    bams = sorted(out_root.glob("bam/*/*.Aligned.sortedByCoord.out.bam"))
    return [Path(b) for b in bams]

def main() -> None:
    """
    Parse CLI args, collect BAM inputs, and run featureCounts end-to-end.

    Either uses user-specified --bam paths or auto-discovers BAMs under out-root.
    Runs featureCounts and writes both raw output and a cleaned PSI-Sigma matrix
    under <out-root>/counts, then prints a short completion summary.
    """
    p = argparse.ArgumentParser(
        description="Run featureCounts and generate a clean gene counts matrix for Psi-Sigma."
    )
    p.add_argument("--threads", type=int, default=12)
    p.add_argument("--stranded", type=int, default=0, choices=[0, 1, 2], help="featureCounts -s (0/1/2)")
    p.add_argument("--paired-end", action="store_true", help="Enable featureCounts paired-end flags: -p -B -C")
    p.add_argument(
        "--bam",
        type=Path,
        action="append",
        default=None,
        help="Optional: provide BAM(s) explicitly. If omitted, auto-discovers under out-root/bam/*/*.Aligned.sortedByCoord.out.bam",
    )
    args = p.parse_args()

    out_dir = OUTDIR
    ensure_dir(out_dir)

    bams = args.bam if args.bam else find_star_bams(out_dir)
    if not bams:
        raise SystemExit(
            f"ERROR: No BAMs found.\n"
            f"Either pass --bam multiple times, or ensure BAMs exist under:\n"
            f"  {out_dir/'bam'}/*/*.Aligned.sortedByCoord.out.bam"
        )

    matrix = featurecounts_gene_counts(
        bams=bams,
        annotation_gtf=ANNOTATION_GTF,
        out_root=out_dir,
        threads=args.threads,
        stranded=args.stranded,
        paired_end=args.paired_end,
    )

    print("DONE")
    print(f"  BAMs counted : {len(bams)}")
    print(f"  Matrix       : {matrix}")


if __name__ == "__main__":
    main()