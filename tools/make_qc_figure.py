#!/usr/bin/env python3
"""
make_qc_figure.py — publication-quality Figure 1 (QC metrics).

WHY NOT THE DEFAULT SCANPY PLOTS
`sc.pl.violin` with jitter draws one point per cell. At 86,000 cells the points
form a solid black mass that completely obscures the violin shape, so the plot
shows nothing about the distribution -- which is the only reason to draw a
violin. It also groups by `condition` (4 groups), hiding the ~2x sequencing
depth variation BETWEEN samples, which matters for interpreting this dataset.

WHAT THIS DOES DIFFERENTLY
  - violins drawn with no jitter, with median and quartile markers
  - grouped by SAMPLE (8), coloured by treatment, so per-sample depth
    heterogeneity is visible
  - total_counts on a log axis (the distribution spans 500 to 140,000)
  - threshold lines drawn where the filters actually cut
  - scatter on log-log axes with the colour bar clipped at the mt threshold,
    so it is not stretched by a handful of dying cells

Run:  python tools/make_qc_figure.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from matplotlib.lines import Line2D

from config import (
    H5AD_RAW, H5AD_QC, FIG_DIR,
    MIN_GENES_PER_CELL, MAX_GENES_PER_CELL, MAX_PCT_MT,
)

sc.settings.verbosity = 0

COL = {"Sucrose": "#4C72B0", "Cocaine": "#C44E52"}
ORDER = ["Female_Sucrose_R1", "Female_Sucrose_R2",
         "Female_Cocaine_R1", "Female_Cocaine_R2",
         "Male_Sucrose_R1", "Male_Sucrose_R2",
         "Male_Cocaine_R1", "Male_Cocaine_R2"]
SHORT = ["F Suc 1", "F Suc 2", "F Coc 1", "F Coc 2",
         "M Suc 1", "M Suc 2", "M Coc 1", "M Coc 2"]


def violin_panel(ax, obs, col, order, log=False, hlines=(), ylabel=""):
    """One violin per sample, no jitter, median and IQR marked."""
    data, colors = [], []
    for s in order:
        v = obs.loc[obs["sample"] == s, col].values
        v = v[v > 0] if log else v
        data.append(np.log10(v) if log else v)
        colors.append(COL[obs.loc[obs["sample"] == s, "treatment"].iloc[0]])

    parts = ax.violinplot(data, positions=range(len(order)), widths=0.8,
                          showextrema=False, showmedians=False)
    for body, c in zip(parts["bodies"], colors):
        body.set_facecolor(c)
        body.set_alpha(0.75)
        body.set_edgecolor("black")
        body.set_linewidth(0.5)

    # Median and interquartile range — the numbers a reader actually wants
    for i, d in enumerate(data):
        q1, med, q3 = np.percentile(d, [25, 50, 75])
        ax.vlines(i, q1, q3, color="black", linewidth=3, zorder=3)
        ax.plot(i, med, "o", color="white", markersize=3.5,
                markeredgecolor="black", markeredgewidth=0.5, zorder=4)

    for y in hlines:
        ax.axhline(np.log10(y) if log else y, ls="--", lw=0.8,
                   color="darkred", zorder=1)

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(SHORT, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    if log:
        lo, hi = ax.get_ylim()
        ticks = np.arange(np.floor(lo), np.ceil(hi) + 1)
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"$10^{{{int(t)}}}$" for t in ticks], fontsize=8)


def main():
    for path, label in [(H5AD_RAW, "before"), (H5AD_QC, "after")]:
        if not path.exists():
            print(f"skipping {label}: {path} not found")
            continue

        adata = sc.read_h5ad(path)
        obs = adata.obs
        if "pct_counts_mt" not in obs.columns:
            adata.var["mt"] = adata.var_names.str.startswith("mt:")
            adata.var["ribo"] = adata.var_names.str.startswith(("RpS", "RpL"))
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"],
                                       percent_top=None, log1p=False, inplace=True)
            obs = adata.obs

        order = [s for s in ORDER if s in set(obs["sample"].astype(str))]
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))

        violin_panel(axes[0], obs, "n_genes_by_counts", order,
                     hlines=[h for h in (MIN_GENES_PER_CELL, MAX_GENES_PER_CELL)
                             if h and label == "before"],
                     ylabel="Genes per cell")
        violin_panel(axes[1], obs, "total_counts", order, log=True,
                     ylabel="UMI counts per cell (log$_{10}$)")
        violin_panel(axes[2], obs, "pct_counts_mt", order,
                     hlines=[MAX_PCT_MT] if label == "before" else [],
                     ylabel="Mitochondrial reads (%)")

        axes[0].set_ylim(bottom=0)   # gene counts cannot be negative
        if label == "before":
            axes[2].set_ylim(0, 15)  # the 82% outliers flatten everything else

        handles = [Line2D([0], [0], marker="s", color="none", markersize=8,
                          markerfacecolor=COL[t], markeredgecolor="black", label=t)
                   for t in ["Sucrose", "Cocaine"]]
        axes[0].legend(handles=handles, frameon=False, fontsize=8, loc="upper left")

        n = adata.n_obs
        fig.suptitle(f"Quality control metrics per sample — {label} filtering "
                     f"({n:,} cells)", fontsize=10, y=1.02)
        fig.tight_layout()
        out = FIG_DIR / f"fig1_qc_{label}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {out}")

        if label == "before":
            # Panel D: the count-vs-gene relationship, log-log.
            fig, ax = plt.subplots(figsize=(5, 4.2))
            x = obs["total_counts"].values
            y = obs["n_genes_by_counts"].values
            c = np.clip(obs["pct_counts_mt"].values, 0, MAX_PCT_MT)
            s = ax.scatter(x, y, c=c, s=1.5, alpha=0.4, cmap="viridis",
                           rasterized=True, linewidths=0)
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.axhline(MIN_GENES_PER_CELL, ls="--", lw=0.8, color="darkred")
            if MAX_GENES_PER_CELL:
                ax.axhline(MAX_GENES_PER_CELL, ls="--", lw=0.8, color="darkred")
            ax.set_xlabel("UMI counts per cell", fontsize=9)
            ax.set_ylabel("Genes per cell", fontsize=9)
            ax.tick_params(labelsize=8)
            ax.spines[["top", "right"]].set_visible(False)
            cb = fig.colorbar(s, ax=ax)
            cb.set_label(f"Mitochondrial % (clipped at {MAX_PCT_MT})", fontsize=8)
            cb.ax.tick_params(labelsize=8)
            ax.set_title("Counts vs. genes detected", fontsize=10)
            fig.tight_layout()
            out = FIG_DIR / "fig1_qc_scatter.png"
            fig.savefig(out, dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"wrote {out}")

        del adata

    print("\nFor the report: use fig1_qc_before.png as Figure 1.")
    print("The per-sample grouping shows the ~2x depth difference between")
    print("samples, and panel 3 shows the male cocaine/sucrose mt% gap.")


if __name__ == "__main__":
    main()
