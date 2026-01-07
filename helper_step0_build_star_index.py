
"""
helper_step0_build_star_index.py

Build a STAR genome index (reference) directory, optionally downloading missing
reference inputs (genome FASTA and annotation GTF) from user-provided URLs.

This script is intended to be the "reference preparation" step that you run once
per (genome FASTA, annotation GTF, read length) combination. It will skip work if
the index output directory already appears to contain a STAR index.

Behavior
--------
1) Verifies that STAR is available on PATH.
2) Ensures genome FASTA and annotation GTF exist:
   - If a file is missing and BOTH a download URL and --download-dir are provided,
     it downloads the file into:
        <download-dir>/<basename-of-requested-file>
   - If a file is missing and no URL/download-dir is provided, exits with an error.
3) Creates the index directory if needed.
4) Skips index generation if the sentinel file "<index-dir>/Genome" exists.
5) Runs STAR genome generation:
       STAR --runMode genomeGenerate
            --genomeDir <index-dir>
            --genomeFastaFiles <genome-fasta>
            --sjdbGTFfile <annotation-gtf>
            --sjdbOverhang <read_length - 1>
            --runThreadN <threads>

Important Notes
---------------
- sjdbOverhang should match sequencing read length (typically read_length - 1).
  If you align reads of a different length than the index was built for, junction
  annotation performance can degrade. Keep separate indices per read length.
- This script does not decompress .gz inputs; provide inputs in the format STAR
  expects in your environment or extend the script to decompress.

Command-line Interface (CLI)
----------------------------
Required:
  --genome-fasta PATH
  --annotation-gtf PATH
  --index-dir PATH
  --read-length INT

Optional:
  --threads INT
  --download-dir PATH
  --genome-url URL
  --annotation-url URL

Examples
--------
Build an index from local files:
  python star_get_index.py \\
    --genome-fasta /refs/mm10.fa \\
    --annotation-gtf /refs/gencode.mm10.gtf \\
    --index-dir /indexes/mm10_rl100 \\
    --read-length 100 \\
    --threads 12

Download missing inputs then build:
  python star_get_index.py \\
    --genome-fasta GRCh38.primary_assembly.genome.fa.gz \\
    --annotation-gtf gencode.v44.annotation.gtf.gz \\
    --index-dir /indexes/GRCh38_rl150 \\
    --read-length 150 \\
    --threads 12 \\
    --download-dir /refs \\
    --genome-url "https://..." \\
    --annotation-url "https://..."
"""

#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import argparse
import urllib.request

from utils import ensure_dir, run, which_or_die

def download_if_missing(url: str, dest: Path) -> None:
    ensure_dir(dest.parent)
    if dest.exists():
        return
    print(f"Downloading {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)


def build_star_index(
    genome_fasta: Path,
    annotation_gtf: Path,
    star_index_dir: Path,
    threads: int,
    read_length: int,
    download_dir: Path | None = None,
    genome_url: str | None = None,
    annotation_url: str | None = None,
) -> None:
    which_or_die("STAR")

    # If missing, optionally download to download_dir/<original filename>
    if not genome_fasta.exists():
        if genome_url and download_dir:
            genome_fasta = download_dir / genome_fasta.name
            download_if_missing(genome_url, genome_fasta)
        else:
            raise SystemExit(f"ERROR: genome FASTA not found: {genome_fasta}")

    if not annotation_gtf.exists():
        if annotation_url and download_dir:
            annotation_gtf = download_dir / annotation_gtf.name
            download_if_missing(annotation_url, annotation_gtf)
        else:
            raise SystemExit(f"ERROR: annotation GTF not found: {annotation_gtf}")

    ensure_dir(star_index_dir)

    sentinel = star_index_dir / "Genome"
    if sentinel.exists():
        print(f"STAR index exists, skipping: {star_index_dir}")
        return

    sjdb_overhang = max(read_length - 1, 1)

    run([
        "STAR",
        "--runThreadN", str(threads),
        "--runMode", "genomeGenerate",
        "--genomeDir", str(star_index_dir),
        "--genomeFastaFiles", str(genome_fasta),
        "--sjdbGTFfile", str(annotation_gtf),
        "--sjdbOverhang", str(sjdb_overhang),
    ])


def main() -> None:
    p = argparse.ArgumentParser(description="Build STAR genome index (optionally download FASTA/GTF).")
    p.add_argument("--genome-fasta", type=Path, required=True)
    p.add_argument("--annotation-gtf", type=Path, required=True)
    p.add_argument("--index-dir", type=Path, required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--read-length", type=int, required=True)

    p.add_argument("--download-dir", type=Path, default=None,
                   help="If provided, missing FASTA/GTF will be downloaded here (requires URLs).")
    p.add_argument("--genome-url", type=str, default=None)
    p.add_argument("--annotation-url", type=str, default=None)

    args = p.parse_args()

    build_star_index(
        genome_fasta=args.genome_fasta,
        annotation_gtf=args.annotation_gtf,
        star_index_dir=args.index_dir,
        threads=args.threads,
        read_length=args.read_length,
        download_dir=args.download_dir,
        genome_url=args.genome_url,
        annotation_url=args.annotation_url,
    )


if __name__ == "__main__":
    main()
