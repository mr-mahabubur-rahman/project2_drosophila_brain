#!/usr/bin/env python3
"""
make_cluster_figure.py — publication-quality Figure 2 (clustering / atlas).

WHY NOT THE DEFAULT SCANPY PLOTS
`sc.pl.umap(color='treatment')` draws all cells in one pass, so whichever
category is plotted LAST covers the others. The result looks as though one
treatment dominates the embedding when the cells are actually interleaved.
The same applies to `condition` and `sample`. Those panels are plotting
artefacts, not data.

The fix is to facet: one panel per category, all cells drawn in grey
underneath so the reader can see where that group sits within the whole.

PANELS
  A  UMAP coloured by Leiden cluster, labelled at cluster centroids
  B  cluster count vs. resolution (the paper's stability argument)
  C  UMAP faceted by sex   -- tests whether sex drives the embedding
  D  UMAP faceted by treatment
  E  per-cluster sample composition, as a stacked bar -- the quantitative
     version of "no cluster is dominated by one sample"

Run:  python tools/make_cluster_figure.py
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

from config import H5AD_CLUSTERED, FIG_DIR, TABLE_DIR, LEIDEN_KEY, TARGET_N_CLUSTERS

sc.settings.verbosity = 0
PT = 0.8          # point size; 86k cells need small points
GREY = "#DDDDDD"


def scatter(ax, xy, mask=None, color=None, cmap_vals=None):
    if mask is not None:
        ax.scatter(xy[:, 0], xy[:, 1], s=PT, c=GREY, linewidths=0, rasterized=True)
        ax.scatter(xy[mask, 0], xy[mask, 1], s=PT, c=color, linewidths=0,
                   rasterized=True)
    else:
        ax.scatter(xy[:, 0], xy[:, 1], s=PT, c=cmap_vals, linewidths=0,
                   rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def main():
    adata = sc.read_h5ad(H5AD_CLUSTERED)
    xy = adata.obsm["X_umap"]
    clusters = adata.obs[LEIDEN_KEY].astype(str).values
    uniq = sorted(set(clusters), key=int)
    print(f"{adata.n_obs:,} cells, {len(uniq)} clusters")

    # ---- A: clusters, labelled at centroids -----------------------------
    fig, ax = plt.subplots(figsize=(7.5, 7))
    cmap = plt.get_cmap("tab20")
    for i, cl in enumerate(uniq):
        m = clusters == cl
        ax.scatter(xy[m, 0], xy[m, 1], s=PT, color=cmap(i % 20),
                   linewidths=0, rasterized=True)
    for cl in uniq:
        m = clusters == cl
        cx, cy = np.median(xy[m, 0]), np.median(xy[m, 1])
        ax.text(cx, cy, cl, fontsize=9, fontweight="bold", ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none",
                          alpha=0.75))
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f"Leiden clustering, resolution 0.8 — {len(uniq)} clusters",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2a_umap_clusters.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote fig2a_umap_clusters.png")

    # ---- C/D: faceted metadata ------------------------------------------
    # One panel per category, rest of the cells in grey. This is the panel
    # that the default overlay plot gets wrong.
    for key, fname in [("sex", "fig2c_umap_by_sex"),
                       ("treatment", "fig2d_umap_by_treatment"),
                       ("sample", "fig2e_umap_by_sample")]:
        cats = list(adata.obs[key].cat.categories)
        ncol = min(4, len(cats))
        nrow = int(np.ceil(len(cats) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3 * ncol, 3 * nrow))
        axes = np.atleast_1d(axes).ravel()
        vals = adata.obs[key].astype(str).values
        for i, c in enumerate(cats):
            m = vals == c
            scatter(axes[i], xy, mask=m, color="#C44E52")
            axes[i].set_title(f"{c}\n({m.sum():,} cells)", fontsize=9)
        for j in range(len(cats), len(axes)):
            axes[j].axis("off")
        fig.suptitle(f"UMAP by {key} — each panel highlights one group",
                     fontsize=11, y=1.01)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"{fname}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {fname}.png")

    # ---- B: resolution sweep ---------------------------------------------
    sweep_path = TABLE_DIR / "resolution_sweep.csv"
    if sweep_path.exists():
        sw = pd.read_csv(sweep_path)
        fig, ax = plt.subplots(figsize=(4.5, 3.4))
        ax.plot(sw.resolution, sw.n_clusters, "o-", color="#4C72B0", lw=1.8, ms=6)
        ax.axhline(TARGET_N_CLUSTERS, ls="--", color="#C44E52", lw=1.2,
                   label=f"Baker et al.: {TARGET_N_CLUSTERS}")
        ax.axvline(0.8, ls=":", color="grey", lw=1)
        ax.annotate("resolution used\nby the paper", (0.8, sw.n_clusters.min()),
                    xytext=(0.85, sw.n_clusters.min() + 1), fontsize=7, color="grey")
        ax.set_xlabel("Leiden resolution", fontsize=9)
        ax.set_ylabel("Number of clusters", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig2b_resolution_sweep.png", dpi=300,
                    bbox_inches="tight")
        plt.close(fig)
        print("wrote fig2b_resolution_sweep.png")

    # ---- E: cluster composition, quantitative ---------------------------
    # The paper asserts no cluster is dominated by one sample. This shows it
    # as numbers rather than as a colour wash.
    comp = pd.crosstab(adata.obs[LEIDEN_KEY], adata.obs["sample"], normalize="index")
    comp = comp.loc[sorted(comp.index, key=lambda x: int(x))]

    fig, ax = plt.subplots(figsize=(10, 4))
    bottom = np.zeros(len(comp))
    cmap2 = plt.get_cmap("tab10")
    for i, s in enumerate(comp.columns):
        ax.bar(range(len(comp)), comp[s].values, bottom=bottom, width=0.85,
               color=cmap2(i), label=s.replace("_", " "))
        bottom += comp[s].values
    ax.axhline(1 / len(comp.columns), ls="--", color="black", lw=1)
    ax.text(len(comp) - 0.4, 1 / len(comp.columns) + 0.015,
            f"even mixing = {100/len(comp.columns):.1f}%", fontsize=7, ha="right")
    ax.set_xticks(range(len(comp)))
    ax.set_xticklabels(comp.index, fontsize=7)
    ax.set_xlim(-0.6, len(comp) - 0.4)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Leiden cluster", fontsize=9)
    ax.set_ylabel("Fraction of cells", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7, ncol=4, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.18))
    ax.set_title("Sample composition per cluster", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2f_cluster_composition.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("wrote fig2f_cluster_composition.png")

    mx = comp.max(axis=1)
    print(f"\nMost sample-skewed cluster: {mx.idxmax()} at {mx.max():.1%} "
          f"(even mixing = {100/len(comp.columns):.1f}%)")
    print(f"Clusters above 30% from one sample: "
          f"{sorted(mx[mx > 0.3].index, key=int)}")


if __name__ == "__main__":
    main()
