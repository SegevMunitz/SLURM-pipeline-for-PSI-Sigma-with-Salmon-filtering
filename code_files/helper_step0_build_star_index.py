#!/usr/bin/env python3
"""
helper_step0_build_star_index.py

Prepare a STAR genome reference index for RNA-seq alignment.

This script builds a STAR genome index from a reference genome FASTA and
an annotation GTF file, optionally downloading missing inputs from user-provided
URLs. It is designed to be run once per reference genome and read length, and
the resulting index can then be reused across multiple alignment jobs.

The script performs the following high-level steps:
  1) Verifies that the STAR executable is available in the environment.
  2) Checks that the genome FASTA and annotation GTF files exist locally.
     - If a required file is missing and both a download URL and a download
       directory are provided, the file is downloaded automatically.
     - If a required file is missing and no download information is provided,
       the script exits with an error.
  3) Creates the STAR index output directory if it does not already exist.
  4) Detects whether a STAR index is already present (via the 'Genome' sentinel
     file) and skips index generation if so.
  5) Runs STAR in genomeGenerate mode using parameters appropriate for the
     specified read length.

The STAR index is generated with a splice junction database overhang
(sjdbOverhang) of (read_length - 1), which should match the read length of the
RNA-seq data that will be aligned. For best results, separate STAR indices
should be generated for datasets with different read lengths.

This script does not modify or decompress input files; genome FASTA and GTF
inputs must be provided in a format supported by STAR in the target environment.

Typical usage is as a preparatory step in an RNA-seq pipeline, executed prior
to any alignment jobs and reused across multiple samples and Slurm array runs.
"""


from __future__ import annotations

from pathlib import Path
import argparse
import urllib.request
import shutil
import gzip

from utils import ensure_dir, run, which_or_die

def download_if_missing(url: str, dest: Path) -> None:
    """
    Download a file from a URL if it does not already exist.

    If URL ends with .gz and dest does NOT end with .gz, download to dest.gz
    and decompress to dest.
    """
    ensure_dir(dest.parent)

    # Case: want uncompressed dest but URL is gz
    if url.endswith(".gz") and not str(dest).endswith(".gz"):
        gz_path = dest.with_suffix(dest.suffix + ".gz")  # e.g. .fa -> .fa.gz
        if dest.exists():
            return
        if not gz_path.exists():
            print(f"Downloading {url} -> {gz_path}")
            urllib.request.urlretrieve(url, gz_path)

        print(f"Decompressing {gz_path} -> {dest}")
        with gzip.open(gz_path, "rb") as fin, dest.open("wb") as fout:
            shutil.copyfileobj(fin, fout)
        return

    # Normal download (dest matches URL)
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
    """
    Build a STAR genome index for a given genome and annotation.

    Verifies STAR availability, ensures FASTA and GTF inputs exist (optionally
    downloading them), and runs STAR genomeGenerate with parameters appropriate
    for the specified read length. Skips index creation if an existing index
    is detected in the output directory.
    """
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
    """
    Parse command-line arguments and invoke STAR index generation.
    
    Handles user input validation and forwards all parameters to
    build_star_index(), which performs the actual reference preparation.
    Intended to be called once per genome / annotation / read-length setup.
    """
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
