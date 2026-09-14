#!/usr/bin/env python3
"""
make_annotation_figure.py — publication-quality Figure 3 (cell-type annotation).

TWO PROBLEMS WITH THE DEFAULTS

1. `sc.pl.umap(color='cell_type', legend_loc='on data')` writes every label at
   its cluster centroid. With 30 clusters and long names ("Central brain
   GABAergic neurons type 1") the labels overlap into an unreadable pile.
   Fix: number the clusters on the plot, put the names in a side legend.

2. The dotplot labels rows by cluster NUMBER, which tells a reader nothing, and
   shows gene names as they appear in this dataset's reference -- which
   substitutes vertebrate ortholog names for five fly symbols. A Drosophila
   figure captioned "VGlut1" implies the mammalian paralog. Fix: relabel rows
   as "12 — Kenyon cells" and display fly symbols with the reference name in
   parentheses.

Run:  python tools/make_annotation_figure.py
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
from matplotlib.lines import Line2D

from config import (
    H5AD_ANNOTATED, FIG_DIR, LEIDEN_KEY,
    CANONICAL_MARKERS, REFERENCE_SYMBOL_ALIASES,
)

sc.settings.verbosity = 0


def display_name(g):
    """VGlut1 -> 'VGlut (VGlut1)' so both the fly symbol and what is in the
    data are visible. Unaliased genes are shown unchanged."""
    fly = REFERENCE_SYMBOL_ALIASES.get(g)
    return f"{fly} ({g})" if fly else g


def umap_with_legend(adata):
    xy = adata.obsm["X_umap"]
    ct = adata.obs["cell_type"].astype(str).values
    cl = adata.obs[LEIDEN_KEY].astype(str).values

    # One entry per cluster, ordered by size, so the legend reads sensibly.
    pairs = (pd.DataFrame({"cluster": cl, "cell_type": ct})
             .groupby(["cluster", "cell_type"], observed=True).size()
             .reset_index(name="n").sort_values("n", ascending=False))

    cmap = plt.get_cmap("tab20")
    colors = {r.cluster: cmap(i % 20) for i, r in enumerate(pairs.itertuples())}

    fig, ax = plt.subplots(figsize=(13, 8))
    for c in pairs["cluster"]:
        m = cl == c
        ax.scatter(xy[m, 0], xy[m, 1], s=0.8, color=colors[c], linewidths=0,
                   rasterized=True)
    for c in pairs["cluster"]:
        m = cl == c
        ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), c, fontsize=8,
                fontweight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                          alpha=0.8))

    handles = [Line2D([0], [0], marker="o", color="none", markersize=7,
                      markerfacecolor=colors[r.cluster],
                      label=f"{r.cluster} — {r.cell_type} ({r.n:,})")
               for r in pairs.itertuples()]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=7.5, frameon=False, ncol=1)

    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("Cell-type annotation of 30 Leiden clusters", fontsize=12)
    fig.tight_layout()
    out = FIG_DIR / "fig3a_umap_celltypes.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out.name}")


def dotplot_labelled(adata):
    # Row label: "12 — Kenyon cells". Ordered so glia group together, then
    # the unannotated clusters sit at the bottom rather than scattered.
    lut = (adata.obs.groupby(LEIDEN_KEY, observed=True)["cell_type"]
           .first().astype(str))
    order = sorted(lut.index, key=lambda c: (
        "Unannotated" in lut[c], lut[c], int(c)))
    labels = {c: f"{c} — {lut[c]}" for c in order}

    adata.obs["_label"] = adata.obs[LEIDEN_KEY].map(labels).astype("category")
    adata.obs["_label"] = adata.obs["_label"].cat.reorder_categories(
        [labels[c] for c in order])

    var_groups = {}
    for cell_type, genes in CANONICAL_MARKERS.items():
        present = [g for g in genes if g in adata.raw.var_names]
        if present:
            var_groups[cell_type] = present

    dp = sc.pl.dotplot(
        adata, var_names=var_groups, groupby="_label", use_raw=True,
        standard_scale="var", show=False, figsize=(13, 9),
        colorbar_title="Scaled mean\nexpression",
        size_title="Fraction of cells\nin cluster (%)",
    )
    # Relabel the gene axis with fly symbols. The reference's ortholog names
    # would otherwise imply mammalian paralogs in a Drosophila figure.
    # scanpy >= 1.10 returns the axes dict directly; older versions return a
    # DotPlot object with .get_axes(). Handle both.
    axes_dict = dp if isinstance(dp, dict) else dp.get_axes()
    axd = axes_dict["mainplot_ax"]
    axd.set_xticklabels([display_name(t.get_text()) for t in axd.get_xticklabels()],
                        rotation=90, fontsize=8)
    axd.tick_params(axis="y", labelsize=8)

    out = FIG_DIR / "fig3b_marker_dotplot.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"wrote {out.name}")

    aliased = [g for gs in var_groups.values() for g in gs
               if g in REFERENCE_SYMBOL_ALIASES]
    if aliased:
        print("  gene axis relabelled: " +
              ", ".join(f"{g}->{REFERENCE_SYMBOL_ALIASES[g]}" for g in aliased))


def main():
    adata = sc.read_h5ad(H5AD_ANNOTATED)
    n_ct = adata.obs["cell_type"].nunique()
    print(f"{adata.n_obs:,} cells, {adata.obs[LEIDEN_KEY].nunique()} clusters, "
          f"{n_ct} distinct labels")
    if n_ct <= 1:
        sys.exit("cell_type is not populated — run scripts/04 with the filled worksheet.")

    umap_with_legend(adata)
    dotplot_labelled(adata)

    n_un = adata.obs["cell_type"].astype(str).str.contains("Unannotated").sum()
    print(f"\nUnannotated: {n_un:,} cells ({100*n_un/adata.n_obs:.1f}%)")
    print("Report these as unannotated rather than guessing -- the paper did")
    print("the same for their C16.")


if __name__ == "__main__":
    main()
