#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import argparse
import sys


def bam_path(out_dir: Path, sample: str) -> Path:
    return out_dir / "bam" / sample / f"{sample}.Aligned.sortedByCoord.out.bam"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--samples", type=Path, required=True, help="samples.txt (one sample ID per line)")
    ap.add_argument("--map", type=Path, required=True, help="TSV: <sample_id>\\t<group_name>")
    args = ap.parse_args()

    out_dir = args.out_dir
    groups_dir = out_dir / "groups"
    groups_dir.mkdir(parents=True, exist_ok=True)

    # read samples
    samples = [ln.strip() for ln in args.samples.read_text().splitlines() if ln.strip()]
    sample_set = set(samples)

    # read mapping
    mapping_lines = [ln.strip() for ln in args.map.read_text().splitlines() if ln.strip()]
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
