#!/usr/bin/env python3
"""
step7_psi_sigma_filtered_by_salmon.py

Filter PSI-Sigma results based on Salmon TPM expression levels.
This script processes PSI-Sigma output files for multiple comparisons,
adding mean TPM values for each group and filtering events based on expression
thresholds and statistical criteria.
"""
import pandas as pd
import re
from pathlib import Path
from typing import Union
from paths import (
    OUTDIR,
    SALMON_TPM_THRESHOLD,
    SALMON_FILTER_MODE,
    COMPARISONS,
    PSISIGMA_ABSPSI_MIN,
    PSISIGMA_P_MAX,
    PSISIGMA_FDR_MAX
)

def event_tpm(ref_tx: str, tpm_series: pd.Series) -> float:
    if pd.isna(ref_tx):
        return 0.0
    toks = [strip_version(t) for t in re.split(r"[;,|]\s*", str(ref_tx)) if t]
    vals = [tpm_series.get(t, 0.0) for t in toks]
    return float(max(vals)) if vals else 0.0

def strip_version(x) -> str:
    """Remove version suffix from transcript IDs (e.g., ENST00000000.1 -> ENST00000000)"""
    return str(x).split(".")[0]

def extract_sample_id_from_bam_path(bam_path: Union[str, Path]) -> str:
    """
    Extract sample ID from BAM path.
    
    Expected format: /path_to_bam/bam/SampleID/SampleID.Aligned.sortedByCoord.out.bam
    Returns: SampleID
    """
    path = Path(bam_path)
    # The sample ID is the parent directory name (the directory containing the BAM)
    return path.parent.name

def load_salmon_tpms(sample_ids):
    """Load Salmon TPM values for multiple samples."""
    dfs = []
    for sid in sample_ids:
        qsf = OUTDIR / "salmon" / sid / "quant.sf"
        if not qsf.exists():
            raise FileNotFoundError(
                f"Salmon quantification file not found for sample '{sid}'. "
                f"Expected: {qsf}\n"
                f"Did step6 (Salmon quantification) complete successfully for all samples?"
            )
        
        q = pd.read_csv(qsf, sep="\t", usecols=["Name", "TPM"])
        q["Name"] = q["Name"].map(strip_version)
        q = q.drop_duplicates(subset=["Name"])
        q = q.set_index("Name")["TPM"].rename(sid)
        dfs.append(q)
    return pd.concat(dfs, axis=1)

def mean_tpm(sample_ids):
    """Calculate mean TPM across samples."""
    return load_salmon_tpms(sample_ids).mean(axis=1)

def validate_psisigma_columns(df, filepath):
    """Validate that PSI-Sigma output has all required columns."""
    required_cols = ["Reference_Transcript", "ΔPSI (%)", "T-test p-value", "FDR (BH)"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    
    if missing_cols:
        raise ValueError(
            f"PSI-Sigma output file is missing required columns.\n"
            f"File: {filepath}\n"
            f"Missing columns: {missing_cols}\n"
            f"Available columns: {list(df.columns)}"
        )
    
if __name__ == "__main__":

    print(f"=== Filtering PSI-Sigma results by Salmon TPM ===")
    print(f"TPM threshold: {SALMON_TPM_THRESHOLD}")
    print(f"Filter mode: {SALMON_FILTER_MODE}")
    print()

    for comp in COMPARISONS:
        comp_dir = OUTDIR / comp.name
        if not comp_dir.exists():
            print(f"WARNING: Comparison directory not found: {comp_dir}, skipping {comp.name}.")
            continue

        print(f"Processing: {comp.name}")

        # Find PSI-Sigma output file
        psisigma_files = list(comp_dir.rglob("*.PSIsigma*.txt"))
        if not psisigma_files:
            print(f"  ERROR: No PSI-Sigma output files found in {comp_dir}")
            print(f"  Expected pattern: PSIsigma*.txt")
            continue

        if len(psisigma_files) > 1:
            print(f"  Warning: Multiple PSI-Sigma files found, using: {psisigma_files[0].name}")
        
        ##################### NEED TO FIX THIS ######################
        psisigma_file = max(psisigma_files, key=lambda p: p.stat().st_mtime)
        print(f"  Reading: {psisigma_file.name}")
        
        # Load PSI-Sigma results
        df = pd.read_csv(psisigma_file, sep="\t")
        print(f"  Total junctions: {len(df)}")
        
        # Validate columns
        validate_psisigma_columns(df, psisigma_file)
                
        # Read group files - these contain BAM paths from step3
        # Note: comp.groupa and comp.groupb are defined in paths.py COMPARISONS
        g1_bam_paths = [x.strip() for x in comp.groupa.read_text().splitlines() if x.strip()]
        g2_bam_paths = [x.strip() for x in comp.groupb.read_text().splitlines() if x.strip()]
        
        # Extract sample IDs from BAM paths
        g1_samples = [extract_sample_id_from_bam_path(path) for path in g1_bam_paths]
        g2_samples = [extract_sample_id_from_bam_path(path) for path in g2_bam_paths]
        
        print(f"  Group 1: {len(g1_samples)} samples: {', '.join(g1_samples)}")
        print(f"  Group 2: {len(g2_samples)} samples: {', '.join(g2_samples)}")
        
        # Calculate mean TPM for each group
        try:
            tpm1 = mean_tpm(g1_samples)
            tpm2 = mean_tpm(g2_samples)
        except FileNotFoundError as e:
            print(f"  ERROR: {e}")
            continue

        # Map TPM values to PSI-Sigma transcripts
        df["meanTPM_group1"] = df["Reference_Transcript"].apply(lambda x: event_tpm(x, tpm1))
        df["meanTPM_group2"] = df["Reference_Transcript"].apply(lambda x: event_tpm(x, tpm2)) 

        # Apply TPM filter
        if SALMON_FILTER_MODE == "both":
            mask_tpm = (
                (df["meanTPM_group1"] >= SALMON_TPM_THRESHOLD) &
                (df["meanTPM_group2"] >= SALMON_TPM_THRESHOLD)
            )
            print(f"  Filter: TPM >= {SALMON_TPM_THRESHOLD} in BOTH groups")
        else:  # "either"
            mask_tpm = (
                (df["meanTPM_group1"] >= SALMON_TPM_THRESHOLD) |
                (df["meanTPM_group2"] >= SALMON_TPM_THRESHOLD)
            )
            print(f"  Filter: TPM >= {SALMON_TPM_THRESHOLD} in EITHER group")
        
        # PSI-Sigma statistical filters
        dpsi_col = "ΔPSI (%)"
        p_col    = "T-test p-value"
        fdr_col  = "FDR (BH)"
        
        # Ensure numeric (handle any non-numeric values)
        for c in (dpsi_col, p_col, fdr_col):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        
        mask_stats = (
            (df[dpsi_col].abs() >= PSISIGMA_ABSPSI_MIN) & 
            (df[p_col] <= PSISIGMA_P_MAX) & 
            (df[fdr_col] <= PSISIGMA_FDR_MAX)
        )
        
        print(f"  Junctions passing TPM filter: {mask_tpm.sum()}")
        print(f"  Junctions passing statistical filter: {mask_stats.sum()}")
        
        # Combine filters
        mask = mask_tpm & mask_stats
        df_filt = df.loc[mask].copy()

        print(f"  Junctions passing BOTH filters: {len(df_filt)}")
        
        # Write output files
        out_xlsx = comp_dir / f"{comp.name}.PSI{PSISIGMA_ABSPSI_MIN}_p{PSISIGMA_P_MAX}_FDR{PSISIGMA_FDR_MAX}.salmon_tpm{SALMON_TPM_THRESHOLD}.{SALMON_FILTER_MODE}.xlsx"
        out_tsv  = comp_dir / f"{comp.name}.PSI{PSISIGMA_ABSPSI_MIN}_p{PSISIGMA_P_MAX}_FDR{PSISIGMA_FDR_MAX}.salmon_tpm{SALMON_TPM_THRESHOLD}.{SALMON_FILTER_MODE}.tsv"
        
        df_filt.to_excel(out_xlsx, index=False)
        df_filt.to_csv(out_tsv, sep="\t", index=False)
        
        print(f"  Wrote: {out_xlsx.name}")
        print(f"  Wrote: {out_tsv.name}")
        print()
    
    print("=== DONE ===")
