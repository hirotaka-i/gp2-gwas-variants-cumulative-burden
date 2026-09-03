# Combine PRS and PCs for GP2 dataset (VWB, release 12)
import numpy as np
import pandas as pd

f_workspace = "/home/jupyter/workspace"
project_name = 'prs_r12'
f_out = f"{f_workspace}/nf_files/{project_name}"
scores = ["Leonard_all", "Leonard_noGBA", "Leonard_noLRRK2", "Leonard_noGBA_noLRRK2"]
ancestries = ["AAC", "AFR", "AJ", "AMR", "CAH", "CAS", "EAS", "EUR", "FIN", "MDE", "SAS"]

# ── Section 1: Build GP2 table ────────────────────────────────────────────────
dfs = []
for anc in ancestries:
    base = f"{f_out}/{anc}/temp/step_08_score"

    # sscore: load each score file and merge on FID/IID
    parts = []
    for sc in scores:
        df = pd.read_csv(f"{base}/{sc}.score.sscore", sep="\t").rename(columns={"#FID": "FID"})
        df = df[["FID", "IID", "ALLELE_CT", "DENOM", "SCORE1_AVG", "SCORE1_SUM"]].rename(columns={
            "ALLELE_CT": f"ALLELE_CT_{sc}", "DENOM": f"DENOM_{sc}",
            "SCORE1_AVG": f"SCORE1_AVG_{sc}", "SCORE1_SUM": f"SCORE1_SUM_{sc}",
        })
        parts.append(df)
    df_anc = parts[0]
    for p in parts[1:]:
        df_anc = df_anc.merge(p, on=["FID", "IID"], how="outer")
    df_anc["Ancestry"] = anc

    df_pc = pd.read_csv(
        f"{f_workspace}/gp2_tier2_eu_release12/raw_genotypes/{anc}/{anc}_release12_vwb.eigenvec",
        sep="\t",
    ).rename(columns={"#FID": "FID"})
    df_pc = df_pc[["FID", "IID"] + [f"PC{i}" for i in range(1, 11)]]

    # Union variant columns across all raw files
    df_raw = None
    for sc in scores:
        rp = pd.read_csv(f"{base}/{sc}.score.raw", sep=r"\s+").rename(columns={"#FID": "FID"})
        rp = rp.drop(columns=[c for c in ["PAT", "MAT", "SEX", "PHENOTYPE"] if c in rp.columns])
        if df_raw is None:
            df_raw = rp
        else:
            extra = [c for c in rp.columns if c not in df_raw.columns]
            if extra:
                df_raw = df_raw.merge(rp[["FID", "IID"] + extra], on=["FID", "IID"])

    dfs.append(df_anc.merge(df_pc, on=["FID", "IID"]).merge(df_raw, on=["FID", "IID"]))

df_all = pd.concat(dfs, ignore_index=True)
df_all.to_csv(f"{f_out}/{len(scores)}PRS_PCs.csv", index=False)
print(f"GP2 dataset shape: {df_all.shape}")