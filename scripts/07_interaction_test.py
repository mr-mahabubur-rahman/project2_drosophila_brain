#!/usr/bin/env python3
"""
07_interaction_test.py — formal sex x treatment interaction test.

WHY THIS MATTERS
Sexual dimorphism is currently inferred by running differential expression
separately in each sex and comparing the resulting gene lists: 90 genes in males
versus 14 in females. That is not a test of difference. A gene at p = 0.04 in
males and p = 0.06 in females appears in one list and not the other while being
statistically indistinguishable -- the "difference between significant and
non-significant is not itself significant" fallacy (Gelman & Stern 2006,
The American Statistician 60:328-331).

Baker et al. use the same approach: separate per-sex analyses compared by Venn
diagram. Neither analysis, as published, tests whether any gene responds
DIFFERENTLY to cocaine in males and females.

WHAT THIS DOES
Fits, per gene:

    expression ~ sex + treatment + sex:treatment

and tests the interaction coefficient. A significant interaction means the
cocaine effect genuinely differs by sex -- which is the claim the dimorphism
result rests on.

WHY PSEUDOBULK RATHER THAN CELL LEVEL
Fitting this at cell level would inherit the pseudoreplication problem: ~86,000
cells but 8 biological replicates. With 8 samples and 4 model terms there are 4
residual degrees of freedom, so power is low and few genes will reach
significance. That is an honest reflection of the design, not a failure of the
method. Report the ranking and the effect sizes; treat individual p-values with
appropriate caution.

The same test can be run per cluster (--per-cluster), which uses more of the
data at the cost of many more tests.

Run:  python scripts/07_interaction_test.py
      python scripts/07_interaction_test.py --per-cluster
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from config import H5AD_ANNOTATED, H5AD_QC, TABLE_DIR, FIG_DIR, PADJ_THRESHOLD

sc.settings.verbosity = 0
MIN_TOTAL_COUNTS = 50


def pseudobulk(adata, cells_mask=None, counts_source=None):
    """
    Sum raw counts per sample -> one profile per biological replicate.

    NOTE ON GENE COVERAGE: after step 03 the `counts` layer holds only the 2,000
    highly variable genes, so using it would silently restrict the test to those.
    `counts_source` supplies the full-gene count matrix from the step-02
    checkpoint; pass it to test all 12,189 genes.
    """
    a = adata[cells_mask] if cells_mask is not None else adata
    if counts_source is not None:
        cs = counts_source[a.obs_names]
        counts, genes = cs.X, counts_source.var_names
    else:
        counts = a.layers["counts"] if "counts" in a.layers else a.raw[a.obs_names].X
        genes = a.var_names if "counts" in a.layers else a.raw.var_names

    rows, meta = [], []
    for s in a.obs["sample"].cat.categories:
        m = (a.obs["sample"] == s).values
        if m.sum() < 10:
            continue
        rows.append(np.asarray(counts[m].sum(axis=0)).ravel())
        o = a.obs.loc[m].iloc[0]
        meta.append({"sample": s, "sex": o["sex"], "treatment": o["treatment"],
                     "n_cells": int(m.sum())})
    if len(rows) < 6:
        return None, None
    pb = pd.DataFrame(np.vstack(rows), index=[m["sample"] for m in meta], columns=genes)
    return pb, pd.DataFrame(meta).set_index("sample")


def fit_interaction(pb, meta):
    """Per-gene OLS of log2 CPM on sex + treatment + interaction."""
    pb = pb.loc[:, pb.sum(axis=0) >= MIN_TOTAL_COUNTS]
    lcpm = np.log2(pb.div(pb.sum(axis=1), axis=0) * 1e6 + 1)

    # Effect coding: male = +0.5, female = -0.5; cocaine = +0.5, sucrose = -0.5.
    # With this coding the interaction term is the difference in cocaine effect
    # between sexes, and the main effects are averages rather than
    # level-specific, which is what we want to interpret.
    sex = np.where(meta.loc[lcpm.index, "sex"] == "Male", 0.5, -0.5)
    trt = np.where(meta.loc[lcpm.index, "treatment"] == "Cocaine", 0.5, -0.5)
    X = sm.add_constant(np.column_stack([sex, trt, sex * trt]))

    rows = []
    for g in lcpm.columns:
        y = lcpm[g].values
        if np.allclose(y, y[0]):
            continue
        try:
            fit = sm.OLS(y, X).fit()
        except Exception:
            continue
        rows.append({
            "gene": g,
            "sex_effect": fit.params[1],
            "treatment_effect": fit.params[2],
            "interaction": fit.params[3],
            "interaction_p": fit.pvalues[3],
            "treatment_p": fit.pvalues[2],
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["interaction_padj"] = multipletests(df.interaction_p, method="fdr_bh")[1]
    df["treatment_padj"] = multipletests(df.treatment_p, method="fdr_bh")[1]
    return df.sort_values("interaction_p")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-cluster", action="store_true")
    args = ap.parse_args()

    adata = sc.read_h5ad(H5AD_ANNOTATED)
    print(f"Loaded {adata.n_obs:,} cells")

    # Full-gene counts from step 02 -- the counts layer here is HVG-only.
    counts_source = None
    if H5AD_QC.exists():
        qc = sc.read_h5ad(H5AD_QC)
        counts_source = qc[adata.obs_names].copy()
        print(f"Using full gene set from step 02: {counts_source.n_vars:,} genes")
        del qc
    else:
        print("WARNING: step-02 checkpoint missing; restricted to HVGs only.")
    print()

    # ---- Global ---------------------------------------------------------
    pb, meta = pseudobulk(adata, counts_source=counts_source)
    print("Design:")
    print(pd.crosstab(meta.sex, meta.treatment).to_string())

    res = fit_interaction(pb, meta)
    res.to_csv(TABLE_DIR / "interaction_global.csv", index=False)

    n_int = int((res.interaction_padj < PADJ_THRESHOLD).sum())
    n_trt = int((res.treatment_padj < PADJ_THRESHOLD).sum())
    print(f"\nGenes tested: {len(res):,}")
    print(f"  significant TREATMENT effect:   {n_trt}")
    print(f"  significant INTERACTION effect: {n_int}")
    print("\nWith 8 samples and 4 model terms there are 4 residual df, so few")
    print("genes are expected to reach significance. The ranking is still")
    print("informative; report effect sizes alongside p-values.")

    print("\nTop 15 by interaction p-value:")
    cols = ["gene", "treatment_effect", "interaction", "interaction_p", "interaction_padj"]
    print(res.head(15)[cols].round(4).to_string(index=False))

    # Volcano of the interaction term
    d = res.copy()
    d["y"] = -np.log10(d.interaction_p.clip(lower=1e-300))
    sig = d.interaction_padj < PADJ_THRESHOLD
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(d.interaction[~sig], d.y[~sig], s=5, c="lightgrey", rasterized=True)
    ax.scatter(d.interaction[sig], d.y[sig], s=18, c="crimson",
               label=f"FDR < 0.05 ({int(sig.sum())})")
    for _, r in d.head(10).iterrows():
        ax.annotate(r["gene"], (r["interaction"], r["y"]), fontsize=7,
                    xytext=(4, 2), textcoords="offset points")
    ax.axvline(0, lw=.6, c="black")
    ax.set_xlabel("Interaction coefficient\n(difference in cocaine effect, male − female)",
                  fontsize=9)
    ax.set_ylabel("−log10 interaction p", fontsize=9)
    ax.set_title("Sex × treatment interaction (pseudobulk, n = 8)", fontsize=11)
    if sig.any():
        ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig5_interaction_volcano.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote fig5_interaction_volcano.png")

    # ---- Compare against the separate-analysis gene lists ----------------
    # The key question: do the genes called male-specific by separate testing
    # actually show a significant interaction?
    male_p = TABLE_DIR / "de_global_male_significant.csv"
    if male_p.exists():
        male_sig = set(pd.read_csv(male_p)["names"])
        sub = res[res.gene.isin(male_sig)]
        if len(sub):
            n = int((sub.interaction_padj < PADJ_THRESHOLD).sum())
            print(f"\nOf {len(sub)} genes called significant in males by separate")
            print(f"per-sex testing, {n} show a significant sex x treatment")
            print("interaction at FDR < 0.05.")
            print("A low number means the separate-analysis approach identifies")
            print("genes that respond in males, not genes that respond")
            print("DIFFERENTLY in males -- a weaker claim than dimorphism.")

    # ---- Per cluster -----------------------------------------------------
    if args.per_cluster:
        print("\n" + "=" * 60)
        print("Per-cluster interaction tests")
        print("=" * 60)
        out = []
        for ct in adata.obs["cell_type"].cat.categories:
            if "excluded" in str(ct):
                continue
            m = (adata.obs["cell_type"] == ct).values
            pbc, metac = pseudobulk(adata, m, counts_source=counts_source)
            if pbc is None:
                continue
            r = fit_interaction(pbc, metac)
            if r.empty:
                continue
            r["cell_type"] = ct
            out.append(r)
            print(f"  {str(ct)[:42]:<42} "
                  f"{int((r.interaction_padj < PADJ_THRESHOLD).sum()):>4} interaction genes")
        if out:
            allr = pd.concat(out, ignore_index=True)
            allr.to_csv(TABLE_DIR / "interaction_per_cluster.csv", index=False)
            print(f"\nwrote interaction_per_cluster.csv ({len(allr):,} rows)")


if __name__ == "__main__":
    main()
