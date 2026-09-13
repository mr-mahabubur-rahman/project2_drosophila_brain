"""
05b_pseudobulk_check.py — validate the cell-level DE with a replicate-level test.

WHY THIS EXISTS (and why it will win you marks)
Script 05 ran Wilcoxon across ~80,000 cells. That test assumes cells are
independent replicates. They are not: all cells in Male_Cocaine_1 came from
one pool of 20 brains processed together. Treating them as 10,000 independent
observations inflates significance -- this is pseudoreplication, and it is the
single most-cited criticism of scRNA-seq DE analyses.

The honest cross-check is PSEUDOBULK: sum counts per sample, giving one
expression profile per biological replicate, then test 2 vs 2. With n=2 per
group this is badly underpowered, so we do NOT expect many significant genes.
What we check is whether the DIRECTION of change agrees with script 05. High
concordance means the cell-level result reflects real biology rather than
pseudoreplication artefact.

Being able to say this in your Discussion -- and in your peer-review rebuttal
-- shows you understand the statistics rather than just the API.

Run:  uv run scripts/05b_pseudobulk_check.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

from config import (
    H5AD_ANNOTATED, FIG_DIR, TABLE_DIR, LOG2FC_THRESHOLD, PADJ_THRESHOLD,
)

sc.settings.verbosity = 1


def make_pseudobulk(adata):
    """Sum raw counts within each sample -> one profile per biological replicate."""
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.raw.X
    genes = adata.raw.var_names if "counts" not in adata.layers else adata.var_names

    rows, meta = [], []
    for sample in adata.obs["sample"].cat.categories:
        mask = (adata.obs["sample"] == sample).values
        if mask.sum() == 0:
            continue
        rows.append(np.asarray(counts[mask].sum(axis=0)).ravel())
        o = adata.obs.loc[mask].iloc[0]
        meta.append({"sample": sample, "sex": o["sex"],
                     "treatment": o["treatment"], "n_cells": int(mask.sum())})

    pb = pd.DataFrame(np.vstack(rows), index=[m["sample"] for m in meta], columns=genes)
    return pb, pd.DataFrame(meta).set_index("sample")


def cpm_log(pb):
    """Counts per million, log2. Crude but adequate for a directional check."""
    return np.log2(pb.div(pb.sum(axis=1), axis=0) * 1e6 + 1)


def main():
    adata = sc.read_h5ad(H5AD_ANNOTATED)
    pb, meta = make_pseudobulk(adata)
    print(f"Pseudobulk matrix: {pb.shape[0]} samples x {pb.shape[1]:,} genes")
    print(meta.to_string())

    # Drop genes with negligible total signal -- they only add noise and
    # multiple-testing burden.
    keep = pb.sum(axis=0) >= 50
    pb = pb.loc[:, keep]
    print(f"Kept {pb.shape[1]:,} genes with >=50 total counts")

    lcpm = cpm_log(pb)
    out = []

    for sex in ["Male", "Female"]:
        idx = meta.index[meta["sex"] == sex]
        coc = [s for s in idx if meta.loc[s, "treatment"] == "Cocaine"]
        suc = [s for s in idx if meta.loc[s, "treatment"] == "Sucrose"]
        if len(coc) < 2 or len(suc) < 2:
            print(f"{sex}: not enough replicates, skipping")
            continue

        a, b = lcpm.loc[coc], lcpm.loc[suc]
        lfc = a.mean(axis=0) - b.mean(axis=0)
        t, p = stats.ttest_ind(a.values, b.values, axis=0, equal_var=True)

        # Benjamini-Hochberg
        p = np.nan_to_num(p, nan=1.0)
        order = np.argsort(p)
        ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
        padj = np.empty_like(ranked)
        padj[order] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)

        df = pd.DataFrame({
            "names": lcpm.columns, "log2FC_pseudobulk": lfc.values,
            "pval": p, "pvals_adj": padj, "sex": sex,
        }).sort_values("pval")
        df.to_csv(TABLE_DIR / f"pseudobulk_de_{sex.lower()}.csv", index=False)

        n_sig = int(((df["pvals_adj"] < PADJ_THRESHOLD) &
                     (df["log2FC_pseudobulk"].abs() > LOG2FC_THRESHOLD)).sum())
        print(f"\n{sex}: {n_sig} genes significant at pseudobulk level (n=2 vs 2)")
        print("  A small number here is EXPECTED, not a failure -- n=2 has almost")
        print("  no power. The concordance check below is the informative part.")
        out.append(df)

    # ---- Concordance with the cell-level result --------------------------
    print("\n--- Direction concordance with cell-level Wilcoxon ---")
    for df in out:
        sex = df["sex"].iloc[0]
        cell_path = TABLE_DIR / f"de_global_{sex.lower()}_all.csv"
        if not cell_path.exists():
            print(f"{sex}: run 05_de_analysis.py first.")
            continue
        cell = pd.read_csv(cell_path)
        merged = cell.merge(df, on="names", suffixes=("_cell", "_pb"))

        # Restrict to genes the cell-level test called significant -- asking
        # whether pseudobulk agrees on the genes we actually claim.
        hits = merged[(merged["pvals_adj_cell"] < PADJ_THRESHOLD) &
                      (merged["logfoldchanges"].abs() > LOG2FC_THRESHOLD)]
        if len(hits) < 5:
            print(f"{sex}: too few cell-level hits to assess.")
            continue

        agree = (np.sign(hits["logfoldchanges"]) ==
                 np.sign(hits["log2FC_pseudobulk"])).mean()
        rho, prho = stats.spearmanr(hits["logfoldchanges"], hits["log2FC_pseudobulk"])
        print(f"{sex}: {len(hits)} cell-level hits | direction agreement "
              f"{agree:.1%} | Spearman rho = {rho:.3f} (p = {prho:.2e})")

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(hits["logfoldchanges"], hits["log2FC_pseudobulk"], s=14, alpha=0.7)
        ax.axhline(0, lw=0.5, c="black"); ax.axvline(0, lw=0.5, c="black")
        ax.set_xlabel("log2FC (cell-level Wilcoxon)")
        ax.set_ylabel("log2FC (pseudobulk, n=2v2)")
        ax.set_title(f"{sex}: concordance, rho = {rho:.2f}", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"05b_pseudobulk_concordance_{sex.lower()}.png", dpi=150)
        plt.close(fig)

    print("\nInterpretation guide for your Discussion:")
    print("  agreement > 80% -> cell-level result is directionally trustworthy")
    print("  agreement ~ 50% -> the signal is likely pseudoreplication artefact")


if __name__ == "__main__":
    main()
