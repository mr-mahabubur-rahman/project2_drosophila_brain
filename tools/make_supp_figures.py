#!/usr/bin/env python3
"""
make_supp_figures.py — replacements for two unusable default plots.

1. PCA VARIANCE (elbow)
   `sc.pl.pca_variance_ratio` annotates all 50 components, so the labels
   overlap into an unreadable smear, and it leaves the y-axis unlabelled.
   This is the plot that justifies N_PCS, so it has to be legible.
   Replacement: clean scree plot, cumulative variance on a second axis,
   and the chosen N_PCS marked.

2. SEX / TREATMENT COMPOSITION
   The faceted UMAPs (fig2c, fig2d) do not work: each group is ~50% of the
   data, so the grey "all cells" layer underneath is almost entirely covered
   and both panels look like the same red blob. The claim being made --
   that neither sex nor treatment drives the embedding -- is quantitative,
   so plot it quantitatively: the fraction of each cluster belonging to each
   group, against the 50% expectation.

Run:  python tools/make_supp_figures.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from config import H5AD_CLUSTERED, FIG_DIR, LEIDEN_KEY, N_PCS

sc.settings.verbosity = 0


def pca_elbow(adata):
    vr = adata.uns["pca"]["variance_ratio"]
    n = min(50, len(vr))
    x = np.arange(1, n + 1)

    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    ax.plot(x, vr[:n] * 100, "o-", color="#4C72B0", ms=4, lw=1.4)
    ax.axvline(N_PCS, ls="--", color="#C44E52", lw=1.2)
    ax.text(N_PCS + 1, vr[:n].max() * 100 * 0.75,
            f"N_PCS = {N_PCS}\n({vr[:N_PCS].sum()*100:.1f}% of variance)",
            fontsize=8, color="#C44E52")
    ax.set_xlabel("Principal component", fontsize=9)
    ax.set_ylabel("Variance explained (%)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    ax2 = ax.twinx()
    ax2.plot(x, np.cumsum(vr[:n]) * 100, color="grey", lw=1.2, ls=":")
    ax2.set_ylabel("Cumulative variance (%)", fontsize=9, color="grey")
    ax2.tick_params(labelsize=8, colors="grey")
    ax2.spines[["top"]].set_visible(False)

    ax.set_title("PCA variance explained", fontsize=11)
    fig.tight_layout()
    out = FIG_DIR / "figS_pca_elbow.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out.name}")
    print(f"  first {N_PCS} PCs capture {vr[:N_PCS].sum()*100:.1f}% of variance")


def composition_bars(adata):
    """Per-cluster fraction by sex and by treatment, against the 50% line."""
    clusters = sorted(adata.obs[LEIDEN_KEY].astype(str).unique(), key=int)

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    for ax, key, groups, colors in [
        (axes[0], "sex", ["Female", "Male"], ["#DD8452", "#4C72B0"]),
        (axes[1], "treatment", ["Cocaine", "Sucrose"], ["#C44E52", "#4C72B0"]),
    ]:
        comp = pd.crosstab(adata.obs[LEIDEN_KEY], adata.obs[key], normalize="index")
        comp = comp.loc[clusters]
        overall = (adata.obs[key].value_counts(normalize=True))[groups[0]]

        bottom = np.zeros(len(comp))
        for g, c in zip(groups, colors):
            ax.bar(range(len(comp)), comp[g].values, bottom=bottom, width=0.85,
                   color=c, label=g)
            bottom += comp[g].values

        ax.axhline(overall, ls="--", color="black", lw=1.2)
        ax.text(len(comp) - 0.3, overall + 0.02,
                f"dataset overall = {overall:.1%}", fontsize=7, ha="right")
        ax.set_ylim(0, 1)
        ax.set_ylabel(f"Fraction by {key}", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=8, frameon=False, ncol=2, loc="upper left")

        dev = (comp[groups[0]] - overall).abs()
        worst = dev.idxmax()
        print(f"  {key}: largest deviation is cluster {worst} "
              f"({comp.loc[worst, groups[0]]:.1%} {groups[0]}, "
              f"expected {overall:.1%})")

    axes[1].set_xticks(range(len(clusters)))
    axes[1].set_xticklabels(clusters, fontsize=7)
    axes[1].set_xlim(-0.6, len(clusters) - 0.4)
    axes[1].set_xlabel("Leiden cluster", fontsize=9)
    axes[0].set_title("Cluster composition by sex and treatment", fontsize=11)
    fig.tight_layout()
    out = FIG_DIR / "figS_composition_sex_treatment.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out.name}")


def main():
    adata = sc.read_h5ad(H5AD_CLUSTERED)
    print(f"{adata.n_obs:,} cells\n")
    pca_elbow(adata)
    print()
    composition_bars(adata)
    print("\nThese replace pca_variance_ratio_03_pca_variance.png,")
    print("fig2c_umap_by_sex.png and fig2d_umap_by_treatment.png.")


if __name__ == "__main__":
    main()
