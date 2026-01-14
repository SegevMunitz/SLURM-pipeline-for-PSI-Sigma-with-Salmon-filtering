#!/usr/bin/env python3
"""
helper_step1_build_salmon_index.py

Prepare a Salmon transcriptome reference index for RNA-seq quantification.

This script builds a Salmon index from a transcriptome FASTA, optionally using a
decoy-aware strategy (recommended) if a genome FASTA is provided. It can also
optionally download missing inputs from user-provided URLs.

Run once per reference (and k-mer size, if you change it). The resulting index
can then be reused across multiple Salmon quant jobs.

High-level steps:
  1) Verifies that the 'salmon' executable is available.
  2) Checks that transcript FASTA (and optionally genome FASTA) exist locally.
     - If missing and URL+download_dir provided, downloads automatically.
  3) Creates the Salmon index output directory if needed.
  4) Detects whether a Salmon index is already present and skips if so.
  5) Builds either:
     - Decoy-aware index (recommended) if genome FASTA provided
     - Transcript-only index otherwise
"""

from __future__ import annotations

from pathlib import Path
import argparse
import shutil
from helper_step0_build_star_index import download_if_missing

from utils import ensure_dir, run, which_or_die

def salmon_index_exists(index_dir: Path) -> bool:
    """
    Heuristic check for a Salmon index.

    Salmon index directories typically contain:
      - info.json
      - ctable.bin
      - refseq.bin
      - mphf.bin
    """
    return (index_dir / "info.json").exists() and (index_dir / "ctable.bin").exists()


def _write_decoys_from_genome_fasta(genome_fasta: Path, out_decoys: Path) -> None:
    """
    Write decoy names (sequence headers) from a genome FASTA to a file.

    Salmon's decoy-aware indexing expects a list of decoy sequence names, one per line.
    This function extracts headers from the genome FASTA (lines beginning with '>').
    """
    ensure_dir(out_decoys.parent)
    print(f"Extracting decoy names from genome FASTA: {genome_fasta} -> {out_decoys}")

    # Stream parse to avoid loading the whole FASTA.
    with genome_fasta.open("rt") as fin, out_decoys.open("wt") as fout:
        for line in fin:
            if line.startswith(">"):
                # Header can contain spaces; decoy name is the first token.
                name = line[1:].strip().split()[0]
                if name:
                    fout.write(name + "\n")


def _cat_fastas(transcripts_fa: Path, genome_fa: Path, out_fa: Path) -> None:
    """
    Concatenate transcript FASTA + genome FASTA into a single "gentrome" FASTA.
    """
    ensure_dir(out_fa.parent)
    print(f"Creating gentrome FASTA: {out_fa}")
    with out_fa.open("wb") as w:
        for src in (transcripts_fa, genome_fa):
            with src.open("rb") as r:
                shutil.copyfileobj(r, w)


def build_salmon_index(
    transcripts_fasta: Path,
    salmon_index_dir: Path,
    threads: int,
    kmer_size: int | None = None,
    # Optional decoy-aware inputs:
    genome_fasta: Path | None = None,
    # Optional downloading:
    download_dir: Path | None = None,
    transcripts_url: str | None = None,
    genome_url: str | None = None,
) -> None:
    """
    Build a Salmon index from transcriptome FASTA.

    If genome_fasta is provided, builds a decoy-aware index:
      - Extract decoy sequence names from genome FASTA headers
      - Concatenate transcripts + genome -> gentrome FASTA
      - salmon index -t gentrome.fa --decoys decoys.txt -i indexdir

    Otherwise builds transcript-only index:
      - salmon index -t transcripts.fa -i indexdir

    Skips if an index is already detected.
    """
    which_or_die("salmon")

    # Resolve/download transcript FASTA if missing
    if not transcripts_fasta.exists():
        if transcripts_url and download_dir:
            transcripts_fasta = download_dir / transcripts_fasta.name
            download_if_missing(transcripts_url, transcripts_fasta)
        else:
            raise SystemExit(f"ERROR: transcripts FASTA not found: {transcripts_fasta}")

    # Resolve/download genome FASTA if missing (only if requested)
    if genome_fasta is not None and not genome_fasta.exists():
        if genome_url and download_dir:
            genome_fasta = download_dir / genome_fasta.name
            download_if_missing(genome_url, genome_fasta)
        else:
            raise SystemExit(f"ERROR: genome FASTA not found: {genome_fasta}")

    ensure_dir(salmon_index_dir)

    if salmon_index_exists(salmon_index_dir):
        print(f"Salmon index exists, skipping: {salmon_index_dir}")
        return

    cmd = ["salmon", "index", "-p", str(threads), "-i", str(salmon_index_dir)]

    if kmer_size is not None:
        cmd.extend(["-k", str(kmer_size)])

    # Decoy-aware (recommended if genome provided)
    if genome_fasta is not None:
        aux_dir = salmon_index_dir.parent / f"{salmon_index_dir.name}__build_aux"
        ensure_dir(aux_dir)
        gentrome = aux_dir / "gentrome.fa"
        decoys = aux_dir / "decoys.txt"

        _write_decoys_from_genome_fasta(genome_fasta, decoys)
        _cat_fastas(transcripts_fasta, genome_fasta, gentrome)

        cmd.extend(["-t", str(gentrome), "--decoys", str(decoys)])
        print("Building decoy-aware Salmon index...")
    else:
        cmd.extend(["-t", str(transcripts_fasta)])
        print("Building transcript-only Salmon index...")

    run(cmd)


def main() -> None:
    p = argparse.ArgumentParser(
        description="Build Salmon index (optionally decoy-aware; optionally download transcripts/genome FASTA)."
    )
    p.add_argument("--transcripts-fasta", type=Path, required=True,
                   help="Transcriptome FASTA (e.g., GENCODE transcripts.fa).")
    p.add_argument("--index-dir", type=Path, required=True,
                   help="Output directory for the Salmon index.")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--kmer-size", type=int, default=None,
                   help="Optional k-mer size (Salmon default is usually fine).")

    # Optional decoy-aware build
    p.add_argument("--genome-fasta", type=Path, default=None,
                   help="If provided, build a decoy-aware index using this genome FASTA.")

    # Optional download
    p.add_argument("--download-dir", type=Path, default=None,
                   help="If provided, missing FASTAs will be downloaded here (requires URLs).")
    p.add_argument("--transcripts-url", type=str, default=None)
    p.add_argument("--genome-url", type=str, default=None)

    args = p.parse_args()

    build_salmon_index(
        transcripts_fasta=args.transcripts_fasta,
        salmon_index_dir=args.index_dir,
        threads=args.threads,
        kmer_size=args.kmer_size,
        genome_fasta=args.genome_fasta,
        download_dir=args.download_dir,
        transcripts_url=args.transcripts_url,
        genome_url=args.genome_url,
    )


if __name__ == "__main__":
    main()
