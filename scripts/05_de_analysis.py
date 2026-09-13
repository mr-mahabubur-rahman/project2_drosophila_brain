"""
05_de_analysis.py — cocaine vs. sucrose differential expression.

WHAT THIS STEP IS FOR
Research Questions 2 and 3. Three levels of comparison, each answering
something different:

  A. GLOBAL, per sex       -> "males respond more than females" (RQ3 magnitude)
  B. PER CLUSTER, per sex  -> "which cell types respond?" (RQ2)
  C. OVERLAP between sexes -> "do they respond with the SAME genes?" (RQ3 identity)

The project guide only covers (A). But (A) alone cannot answer RQ2, and the
paper's central claim -- that Kenyon cells and glia respond most -- lives
entirely in (B). Skipping it would leave a research question unanswered.

A STATISTICAL CAVEAT YOU SHOULD STATE IN YOUR REPORT
Wilcoxon across cells treats every cell as an independent replicate. They are
not -- cells from one fly brain are correlated, so p-values are anti-
conservative and gene counts are inflated relative to a pseudobulk analysis
(2 replicates per group). The paper has the same limitation (they used MAST
per cell). Report the direction and relative magnitude confidently; treat the
absolute gene counts as approximate. Script 05b does the pseudobulk check.

Run:  uv run scripts/05_de_analysis.py
"""

import warnings
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from config import (
    H5AD_ANNOTATED, FIG_DIR, TABLE_DIR,
    LEIDEN_KEY, DE_METHOD, LOG2FC_THRESHOLD, PAPER_EQUIVALENT_LOG2FC,
    PADJ_THRESHOLD, MIN_CELLS_PER_GROUP_FOR_DE,
    PAPER_CORE_RESPONSE_GENES, PAPER_DISCUSSION_GENES,
)

warnings.filterwarnings("ignore", category=FutureWarning)
sc.settings.verbosity = 1
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)


def de_cocaine_vs_sucrose(subset):
    """
    Wilcoxon rank-sum, Cocaine vs Sucrose, on log-normalized values.

    reference='Sucrose' makes sucrose the denominator, so a positive
    logfoldchange means UP in cocaine. Get this backwards and every
    biological statement in your Discussion inverts.
    """
    sc.tl.rank_genes_groups(
        subset, groupby="treatment", reference="Sucrose",
        method=DE_METHOD, use_raw=True, pts=True,
    )
    df = sc.get.rank_genes_groups_df(subset, group="Cocaine")
    return df.dropna(subset=["logfoldchanges", "pvals_adj"])


def significant(df, lfc=LOG2FC_THRESHOLD):
    return df[(df["logfoldchanges"].abs() > lfc) & (df["pvals_adj"] < PADJ_THRESHOLD)]


def enough_cells(subset):
    """Both arms need enough cells, or Wilcoxon returns noise."""
    counts = subset.obs["treatment"].value_counts()
    return (
        counts.get("Cocaine", 0) >= MIN_CELLS_PER_GROUP_FOR_DE
        and counts.get("Sucrose", 0) >= MIN_CELLS_PER_GROUP_FOR_DE
    )


def volcano(df, title, path, lfc=LOG2FC_THRESHOLD, n_label=12):
    d = df.copy()
    d["neglog10p"] = -np.log10(d["pvals_adj"].clip(lower=1e-300))
    sig = (d["pvals_adj"] < PADJ_THRESHOLD) & (d["logfoldchanges"].abs() > lfc)

    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(d.loc[~sig, "logfoldchanges"], d.loc[~sig, "neglog10p"],
               s=4, c="lightgrey", rasterized=True)
    up = sig & (d["logfoldchanges"] > 0)
    dn = sig & (d["logfoldchanges"] < 0)
    ax.scatter(d.loc[up, "logfoldchanges"], d.loc[up, "neglog10p"],
               s=6, c="crimson", label=f"up ({int(up.sum())})", rasterized=True)
    ax.scatter(d.loc[dn, "logfoldchanges"], d.loc[dn, "neglog10p"],
               s=6, c="steelblue", label=f"down ({int(dn.sum())})", rasterized=True)

    for _, r in d[sig].nlargest(n_label, "neglog10p").iterrows():
        ax.annotate(r["names"], (r["logfoldchanges"], r["neglog10p"]),
                    fontsize=6, alpha=0.85)

    ax.axhline(-np.log10(PADJ_THRESHOLD), ls="--", lw=0.7, c="black")
    ax.axvline(lfc, ls="--", lw=0.7, c="black")
    ax.axvline(-lfc, ls="--", lw=0.7, c="black")
    ax.set_xlabel("log2 fold change (Cocaine / Sucrose)")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    adata = sc.read_h5ad(H5AD_ANNOTATED)
    label_col = "cell_type" if adata.obs["cell_type"].nunique() > 1 else LEIDEN_KEY
    print(f"Loaded {adata.n_obs:,} cells. Labelling clusters by '{label_col}'.")

    results = {}

    # ==================================================================
    # A. GLOBAL DE PER SEX  (RQ3, magnitude)
    # ==================================================================
    print("\n" + "=" * 66)
    print("A. Global differential expression, by sex")
    print("=" * 66)

    summary_rows = []
    for sex in ["Male", "Female"]:
        sub = adata[adata.obs["sex"] == sex].copy()
        df = de_cocaine_vs_sucrose(sub)
        results[sex] = df

        df.to_csv(TABLE_DIR / f"de_global_{sex.lower()}_all.csv", index=False)
        sig_guide = significant(df, LOG2FC_THRESHOLD)
        sig_paper = significant(df, PAPER_EQUIVALENT_LOG2FC)
        sig_guide.to_csv(TABLE_DIR / f"de_global_{sex.lower()}_significant.csv", index=False)

        print(f"\n{sex}: {sub.n_obs:,} cells")
        print(f"  |log2FC|>{LOG2FC_THRESHOLD} & padj<{PADJ_THRESHOLD}: "
              f"{len(sig_guide)} genes "
              f"({int((sig_guide['logfoldchanges']>0).sum())} up, "
              f"{int((sig_guide['logfoldchanges']<0).sum())} down)")
        print(f"  paper-equivalent |ln FC|>1 (= |log2FC|>{PAPER_EQUIVALENT_LOG2FC:.2f}): "
              f"{len(sig_paper)} genes")

        volcano(df, f"{sex}: cocaine vs. sucrose",
                FIG_DIR / f"05_volcano_{sex.lower()}.png")

        summary_rows.append({
            "sex": sex, "n_cells": sub.n_obs,
            "n_sig_guide_threshold": len(sig_guide),
            "n_up": int((sig_guide["logfoldchanges"] > 0).sum()),
            "n_down": int((sig_guide["logfoldchanges"] < 0).sum()),
            "n_sig_paper_threshold": len(sig_paper),
        })
        del sub

    pd.DataFrame(summary_rows).to_csv(TABLE_DIR / "de_global_summary.csv", index=False)

    m, f = summary_rows[0]["n_sig_guide_threshold"], summary_rows[1]["n_sig_guide_threshold"]
    print(f"\nMale:Female DE gene ratio = {m}:{f}"
          + (f" ({m/f:.2f}x)" if f else ""))
    print("Paper reported 691 male vs 322 female (~2.1x) at their own thresholds.")
    print("The RATIO is the reproducible claim; the absolute counts depend on")
    print("normalization, thresholds, and test choice, all of which differ here.")

    # ==================================================================
    # B. PER-CLUSTER DE PER SEX  (RQ2)
    # ==================================================================
    print("\n" + "=" * 66)
    print("B. Per-cluster differential expression, by sex")
    print("=" * 66)

    per_cluster = []
    per_cluster_genes = []
    clusters = sorted(adata.obs[label_col].unique(), key=str)

    for sex in ["Male", "Female"]:
        sex_data = adata[adata.obs["sex"] == sex]
        for cl in clusters:
            sub = sex_data[sex_data.obs[label_col] == cl].copy()
            if not enough_cells(sub):
                per_cluster.append({
                    "sex": sex, "cluster": cl, "n_cells": sub.n_obs,
                    "n_sig": np.nan, "n_up": np.nan, "n_down": np.nan,
                    "note": f"skipped: <{MIN_CELLS_PER_GROUP_FOR_DE} cells in an arm",
                })
                del sub
                continue
            try:
                df = de_cocaine_vs_sucrose(sub)
            except Exception as e:
                per_cluster.append({"sex": sex, "cluster": cl, "n_cells": sub.n_obs,
                                    "n_sig": np.nan, "n_up": np.nan, "n_down": np.nan,
                                    "note": f"failed: {e}"})
                del sub
                continue

            sig = significant(df)
            per_cluster.append({
                "sex": sex, "cluster": cl, "n_cells": sub.n_obs,
                "n_sig": len(sig),
                "n_up": int((sig["logfoldchanges"] > 0).sum()),
                "n_down": int((sig["logfoldchanges"] < 0).sum()),
                "note": "",
            })
            if len(sig):
                s = sig.copy()
                s["sex"], s["cluster"] = sex, cl
                per_cluster_genes.append(s)
            del sub
        print(f"  {sex}: done ({len(clusters)} clusters)")

    pc = pd.DataFrame(per_cluster)
    pc.to_csv(TABLE_DIR / "de_per_cluster_summary.csv", index=False)

    if per_cluster_genes:
        allg = pd.concat(per_cluster_genes, ignore_index=True)
        allg.to_csv(TABLE_DIR / "de_per_cluster_significant_genes.csv", index=False)

    ranked = (
        pc.dropna(subset=["n_sig"])
        .groupby("cluster", observed=True)["n_sig"].sum()
        .sort_values(ascending=False)
    )
    print("\nTop 10 most cocaine-responsive clusters (DE genes summed across sexes):")
    print(ranked.head(10).to_string())
    print("\nCompare against the paper: Kenyon cells (their C11), astrocytes (C17),")
    print("surface glia (C22), and unannotated C16 were their strongest responders.")
    print("Match by CELL TYPE, not by cluster number.")

    # Heatmap of DE burden per cluster per sex -- report Figure 4.
    pivot = pc.pivot_table(index="cluster", columns="sex", values="n_sig").fillna(0)
    fig, ax = plt.subplots(figsize=(5, max(4, 0.25 * len(pivot))))
    im = ax.imshow(pivot.values, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=6)
    ax.set_title("Significant DE genes per cluster", fontsize=10)
    fig.colorbar(im, ax=ax, label="n DE genes")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "05_de_burden_heatmap.png", dpi=150)
    plt.close(fig)

    # ==================================================================
    # C. SEXUAL DIMORPHISM OVERLAP  (RQ3, identity)
    # ==================================================================
    print("\n" + "=" * 66)
    print("C. Overlap between male and female responses")
    print("=" * 66)

    male_sig = significant(results["Male"])
    female_sig = significant(results["Female"])
    ms, fs = set(male_sig["names"]), set(female_sig["names"])
    shared = ms & fs

    print(f"  Male only:   {len(ms - fs)}")
    print(f"  Female only: {len(fs - ms)}")
    print(f"  Shared:      {len(shared)}")

    # The paper's subtler point: some shared genes move in OPPOSITE directions
    # between sexes. That is stronger evidence of dimorphism than counts alone.
    if shared:
        mlfc = male_sig.set_index("names")["logfoldchanges"]
        flfc = female_sig.set_index("names")["logfoldchanges"]
        cmp = pd.DataFrame({
            "male_log2FC": mlfc[list(shared)],
            "female_log2FC": flfc[list(shared)],
        })
        cmp["opposite_direction"] = np.sign(cmp["male_log2FC"]) != np.sign(cmp["female_log2FC"])
        cmp.sort_values("opposite_direction", ascending=False).to_csv(
            TABLE_DIR / "de_shared_genes_direction.csv")
        n_opp = int(cmp["opposite_direction"].sum())
        print(f"  Of the shared genes, {n_opp} change in OPPOSITE directions "
              "between sexes.")

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(cmp["male_log2FC"], cmp["female_log2FC"], s=14,
                   c=np.where(cmp["opposite_direction"], "crimson", "grey"))
        lim = np.nanmax(np.abs(cmp.values[:, :2].astype(float))) * 1.1
        ax.plot([-lim, lim], [-lim, lim], ls="--", lw=0.7, c="black")
        ax.axhline(0, lw=0.5, c="black"); ax.axvline(0, lw=0.5, c="black")
        ax.set_xlabel("male log2FC"); ax.set_ylabel("female log2FC")
        ax.set_title("Shared DE genes: concordance between sexes", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "05_sex_concordance.png", dpi=150)
        plt.close(fig)

    pd.DataFrame({
        "category": ["male_only", "female_only", "shared"],
        "n_genes": [len(ms - fs), len(fs - ms), len(shared)],
        "genes": [", ".join(sorted(ms - fs)), ", ".join(sorted(fs - ms)),
                  ", ".join(sorted(shared))],
    }).to_csv(TABLE_DIR / "de_sex_overlap.csv", index=False)

    # ==================================================================
    # D. POSITIVE CONTROL against the paper's named genes
    # ==================================================================
    print("\n" + "=" * 66)
    print("D. Sanity check: genes the paper names explicitly")
    print("=" * 66)
    print("If none of these move at all, suspect a pipeline error before you")
    print("conclude the paper was wrong.\n")

    watch = (PAPER_CORE_RESPONSE_GENES["up"] + PAPER_CORE_RESPONSE_GENES["down"]
             + PAPER_DISCUSSION_GENES)
    rows = []
    for sex in ["Male", "Female"]:
        d = results[sex].set_index("names")
        for g in watch:
            if g in d.index:
                r = d.loc[g]
                rows.append({
                    "gene": g, "sex": sex,
                    "log2FC": round(float(r["logfoldchanges"]), 3),
                    "padj": float(r["pvals_adj"]),
                    "significant": bool(abs(r["logfoldchanges"]) > LOG2FC_THRESHOLD
                                        and r["pvals_adj"] < PADJ_THRESHOLD),
                    "paper_expectation": (
                        "up" if g in PAPER_CORE_RESPONSE_GENES["up"]
                        else "down" if g in PAPER_CORE_RESPONSE_GENES["down"]
                        else "discussed"
                    ),
                })
    if rows:
        chk = pd.DataFrame(rows)
        chk.to_csv(TABLE_DIR / "paper_gene_check.csv", index=False)
        print(chk.pivot_table(index="gene", columns="sex",
                              values="log2FC").round(2).to_string())

    print(f"\nAll tables written to {TABLE_DIR}")
    print(f"All figures written to {FIG_DIR}")
    plt.close("all")


if __name__ == "__main__":
    main()
