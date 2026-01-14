#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from paths import (
    OUTDIR,
    SALMON_TPM_THRESHOLD,
    SALMON_FILTER_MODE,
)

def strip_version(x):
    return str(x).split(".")[0]

def load_salmon_tpms(sample_ids):
    dfs = []
    for sid in sample_ids:
        qsf = OUTDIR / "salmon" / sid / "quant.sf"
        q = pd.read_csv(qsf, sep="\t", usecols=["Name", "TPM"])
        q["Name"] = q["Name"].map(strip_version)
        q = q.drop_duplicates(subset=["Name"])
        q = q.set_index("Name")["TPM"].rename(sid)
        dfs.append(q)
    return pd.concat(dfs, axis=1)

def mean_tpm(sample_ids):
    return load_salmon_tpms(sample_ids).mean(axis=1)

if __name__ == "__main__":

    # ---- ONLY THIS COMPARISON ----
    comp = Path(
        "/sci/labs/zvika.granot/segev.munitz/"
        "neutrophil_RNAseq/psi_sigma/IL4_vs_LPSIFNG"
    )

    psisigma_file = next(comp.glob("PSIsigma*.txt"))
    df = pd.read_csv(psisigma_file, sep="\t")

    if "Reference_Transcript" not in df.columns:
        raise ValueError("Missing column: Reference_Transcript")

    df["Reference_Transcript"] = df["Reference_Transcript"].map(strip_version)

    g1 = [x.strip() for x in (comp / "group1.txt").read_text().splitlines() if x.strip()]
    g2 = [x.strip() for x in (comp / "group2.txt").read_text().splitlines() if x.strip()]

    tpm1 = mean_tpm(g1)
    tpm2 = mean_tpm(g2)

    df["meanTPM_group1"] = df["Reference_Transcript"].map(tpm1).fillna(0.0)
    df["meanTPM_group2"] = df["Reference_Transcript"].map(tpm2).fillna(0.0)

    if SALMON_FILTER_MODE == "both":
        mask = (
            (df["meanTPM_group1"] >= SALMON_TPM_THRESHOLD) &
            (df["meanTPM_group2"] >= SALMON_TPM_THRESHOLD)
        )
    else:  # "either"
        mask = (
            (df["meanTPM_group1"] >= SALMON_TPM_THRESHOLD) |
            (df["meanTPM_group2"] >= SALMON_TPM_THRESHOLD)
        )

    df_filt = df.loc[mask].copy()

    out_xlsx = comp / f"IL4_vs_LPSIFNG.salmon_tpm{SALMON_TPM_THRESHOLD}.{SALMON_FILTER_MODE}.xlsx"
    out_tsv  = comp / f"IL4_vs_LPSIFNG.salmon_tpm{SALMON_TPM_THRESHOLD}.{SALMON_FILTER_MODE}.tsv"

    df_filt.to_excel(out_xlsx, index=False)
    df_filt.to_csv(out_tsv, sep="\t", index=False)

    print(
        f"IL4_vs_LPSIFNG (N2 vs N1): "
        f"kept {len(df_filt)} / {len(df)} events "
        f"(TPM ≥ {SALMON_TPM_THRESHOLD}, mode={SALMON_FILTER_MODE})"
    )
