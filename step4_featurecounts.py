
"""
step4_featurecounts.py

Gene-level read counting from aligned RNA-seq BAM files using Subread
`featureCounts`, plus conversion of the output into a clean count matrix TSV
suitable for downstream tools (e.g., Psi-Sigma).

This module assumes we have:
  - One BAM file per sample (typically STAR-aligned, coordinate-sorted)
  - A matching gene annotation GTF file for the same genome build

Overview
--------
The main workflow implemented here is:

  BAM files (per sample)
    → featureCounts (gene-level exon counting)
    → featureCounts.gene_counts.txt (raw featureCounts output)
    → gene_counts.matrix.tsv (clean matrix: gene_id + sample columns)

The returned output is the cleaned matrix file path.

Key behaviors
-------------
- Uses featureCounts to count reads overlapping exons (`-t exon`) and summarize
  counts by gene ID (`-g gene_id`).
- Supports strandedness via the `stranded` argument (passed to `-s`).
- Supports paired-end counting when `paired_end=True` by enabling:
    -p  : count fragments instead of reads
    -B  : require both ends to be successfully aligned
    -C  : exclude chimeric fragments

Outputs
-------
All outputs are written under:
  <out_root>/counts/

Files created:
  1) featureCounts.gene_counts.txt
     The standard featureCounts output file, including annotation columns
     (Chr/Start/End/Strand/Length) plus one column per BAM/sample.

  2) gene_counts.matrix.tsv
     A simplified tab-delimited matrix intended for downstream differential
     expression / splicing tools. It contains:
       - a first column named "gene_id"
       - one column per sample (derived from the BAM paths in the featureCounts run)
     and excludes featureCounts annotation columns.

Functions
---------
featurecounts_gene_counts(...)
    Runs featureCounts on a list of BAMs and writes:
      - featureCounts.gene_counts.txt (raw)
      - gene_counts.matrix.tsv (clean)
    Returns the path to the cleaned matrix.

make_clean_counts_matrix(featurecounts_out, clean_out)
    Parses the raw featureCounts output and writes a cleaned matrix:
      - skips initial comment lines beginning with '#'
      - keeps the header row and data rows
      - drops featureCounts annotation columns:
          Chr, Start, End, Strand, Length
      - renames "Geneid" to "gene_id"
    Writes the result as a tab-separated TSV.

Arguments and conventions
-------------------------
- bams: list of BAM file paths to include as columns in the count matrix.
- annotation_gtf: GTF file used for read assignment; must match genome build.
- threads: number of CPU threads for featureCounts (-T).
- stranded: featureCounts strandedness flag (-s):
    0 = unstranded
    1 = stranded
    2 = reversely stranded
- paired_end: if True, enable fragment counting with -p -B -C.

Example usage
-------------
Count genes from STAR BAMs (paired-end, unstranded):
  matrix = featurecounts_gene_counts(
      bams=[Path("sample1.bam"), Path("sample2.bam")],
      annotation_gtf=Path("gencode.gtf"),
      out_root=Path("analysis_out"),
      threads=12,
      stranded=0,
      paired_end=True,
  )
  print("Counts matrix:", matrix)

Notes / assumptions
-------------------
- featureCounts must be installed and discoverable on PATH.
- BAMs should be coordinate-sorted and indexed for best performance (indexing
  is not strictly required by featureCounts but is standard in pipelines).
- The cleaned matrix keeps only gene_id and sample count columns; any additional
  metadata columns from featureCounts are removed to simplify downstream tools.
- Column names for samples come from the featureCounts header (often BAM paths).
  If your downstream tool requires short sample IDs, consider renaming columns
  after generation or extend `make_clean_counts_matrix` to map names.
"""

from __future__ import annotations
import argparse
import csv
from pathlib import Path
from typing import List

from utils import ensure_dir, run, which_or_die

def featurecounts_gene_counts(
    bams: List[Path],
    annotation_gtf: Path,
    out_root: Path,
    threads: int,
    stranded: int,
    paired_end: bool,
) -> Path:
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
    bams = sorted(out_root.glob("bam/*/*.Aligned.sortedByCoord.out.bam"))
    return [Path(b) for b in bams]

def main() -> None:
    p = argparse.ArgumentParser(
        description="Run featureCounts and generate a clean gene counts matrix for Psi-Sigma."
    )
    p.add_argument("--out-root", type=Path, required=True, help="Pipeline output root containing bam/ directory.")
    p.add_argument("--annotation-gtf", type=Path, required=True, help="Annotation GTF used for counting.")
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

    out_root = args.out_root.resolve()
    ensure_dir(out_root)

    bams = args.bam if args.bam else find_star_bams(out_root)
    if not bams:
        raise SystemExit(
            f"ERROR: No BAMs found.\n"
            f"Either pass --bam multiple times, or ensure BAMs exist under:\n"
            f"  {out_root/'bam'}/*/*.Aligned.sortedByCoord.out.bam"
        )

    matrix = featurecounts_gene_counts(
        bams=bams,
        annotation_gtf=args.annotation_gtf,
        out_root=out_root,
        threads=args.threads,
        stranded=args.stranded,
        paired_end=args.paired_end,
    )

    print("DONE")
    print(f"  BAMs counted : {len(bams)}")
    print(f"  Matrix       : {matrix}")


if __name__ == "__main__":
    main()