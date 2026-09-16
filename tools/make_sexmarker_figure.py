#!/usr/bin/env python3
"""
make_sexmarker_figure.py — the labelling discrepancy in one figure.

WHAT IT SHOWS
Panel A: per-sample expression of three male-specific and three female-specific
         markers, with samples ordered by their declared sex. If the labels were
         correct, the male markers would be high in the four male-labelled
         samples. They are not: the split follows treatment.
Panel B: the same data as a two-marker scatter. Each point is one sample,
         coloured by declared treatment and shaped by declared sex. Samples
         separate by colour, not by shape.
Panel C: log2 fold changes for the six sex markers in the two nominal treatment
         contrasts, showing that both contrasts return them among their most
         significant genes.

Panel C reads from the differential expression tables, so run scripts/05 first.

Run:  python tools/make_sexmarker_figure.py
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
from matplotlib.patches import Patch

from config import H5AD_QC, FIG_DIR, TABLE_DIR

sc.settings.verbosity = 0

RED, BLUE, DARK, GREY = "#C44E52", "#4C72B0", "#2B2B2B", "#8C8C8C"

MALE_MARKERS = ["lncRNA:roX1", "lncRNA:roX2", "roX1", "roX2"]
FEMALE_MARKERS = ["Yp1", "Yp2", "Yp3", "Sxl"]

ORDER = ["Female_Sucrose_R1", "Female_Sucrose_R2",
         "Female_Cocaine_R1", "Female_Cocaine_R2",
         "Male_Sucrose_R1", "Male_Sucrose_R2",
         "Male_Cocaine_R1", "Male_Cocaine_R2"]
SHORT = ["F Suc 1", "F Suc 2", "F Coc 1", "F Coc 2",
         "M Suc 1", "M Suc 2", "M Coc 1", "M Coc 2"]


def sample_means(adata, genes):
    d = sc.get.obs_df(adata, keys=genes)
    d["sample"] = adata.obs["sample"].astype(str).values
    return d.groupby("sample", observed=True).mean()


def main():
    adata = sc.read_h5ad(H5AD_QC)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    male_g = [g for g in MALE_MARKERS if g in adata.var_names][:2]
    female_g = [g for g in FEMALE_MARKERS if g in adata.var_names][:3]
    genes = male_g + female_g
    if not male_g or not female_g:
        sys.exit(f"markers not found; available roX/Yp: "
                 f"{[v for v in adata.var_names if 'roX' in v or v.startswith('Yp')]}")

    m = sample_means(adata, genes)
    order = [s for s in ORDER if s in m.index]
    short = [SHORT[ORDER.index(s)] for s in order]
    m = m.loc[order]
    trt = ["Cocaine" if "Cocaine" in s else "Sucrose" for s in order]
    colors = [RED if t == "Cocaine" else BLUE for t in trt]

    fig = plt.figure(figsize=(13, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.55, 0.85, 1.05], wspace=0.30)

    # ---------------- Panel A: grouped bars ----------------
    axA = fig.add_subplot(gs[0])
    n_g = len(genes)
    x = np.arange(len(order))
    w = 0.8 / n_g
    for i, g in enumerate(genes):
        axA.bar(x + i * w - 0.4 + w / 2, m[g].values, width=w * 0.92,
                color=colors, edgecolor="none",
                alpha=1.0 if g in male_g else 0.55)
    for i, g in enumerate(genes):
        axA.text(x[0] + i * w - 0.4 + w / 2, m[g].max() * 1.02 + 0.12,
                 g.replace("lncRNA:", ""), rotation=90, fontsize=6.5,
                 ha="center", va="bottom", color=GREY)
    axA.axvline(3.5, color=DARK, lw=1.0, ls="--")
    axA.text(1.5, axA.get_ylim()[1] * 0.97, "declared female", ha="center",
             fontsize=8.5, color=DARK)
    axA.text(5.5, axA.get_ylim()[1] * 0.97, "declared male", ha="center",
             fontsize=8.5, color=DARK)
    axA.set_xticks(x)
    axA.set_xticklabels(short, rotation=45, ha="right", fontsize=8)
    axA.set_ylabel("Mean log-normalised expression", fontsize=9)
    axA.set_title("A  Sex markers by sample", fontsize=10.5, loc="left")
    axA.spines[["top", "right"]].set_visible(False)
    axA.tick_params(labelsize=8)
    axA.legend(handles=[
        Patch(facecolor=BLUE, label="declared sucrose"),
        Patch(facecolor=RED, label="declared cocaine"),
        Patch(facecolor=GREY, alpha=1.0, label="male markers (solid)"),
        Patch(facecolor=GREY, alpha=0.55, label="female markers (faded)"),
    ], fontsize=7, frameon=False, loc="upper center", ncol=2,
        bbox_to_anchor=(0.5, -0.30))

    # ---------------- Panel B: two-marker scatter ----------------
    axB = fig.add_subplot(gs[1])
    xg, yg = male_g[0], female_g[0]
    for s, c, t in zip(order, colors, trt):
        sex = "Female" if s.startswith("Female") else "Male"
        axB.scatter(m.loc[s, xg], m.loc[s, yg], s=120, c=c,
                    marker="o" if sex == "Female" else "^",
                    edgecolor=DARK, linewidth=0.7, zorder=3)
    axB.set_xlabel(f"{xg.replace('lncRNA:', '')}  (male-specific)", fontsize=8.5)
    axB.set_ylabel(f"{yg}  (female-specific)", fontsize=8.5)
    axB.set_title("B  Samples separate by treatment,\n     not by declared sex",
                  fontsize=10.5, loc="left")
    axB.spines[["top", "right"]].set_visible(False)
    axB.tick_params(labelsize=8)
    axB.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GREY,
               markeredgecolor=DARK, markersize=8, label="declared female"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=GREY,
               markeredgecolor=DARK, markersize=8, label="declared male"),
    ], fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22))

    # ---------------- Panel C: fold changes in both contrasts ----------------
    axC = fig.add_subplot(gs[2])
    rows = []
    for sex, lab in [("male", "male-labelled"), ("female", "female-labelled")]:
        p = TABLE_DIR / f"de_global_{sex}_all.csv"
        if not p.exists():
            continue
        d = pd.read_csv(p).set_index("names")
        for g in genes:
            if g in d.index:
                rows.append({"gene": g.replace("lncRNA:", ""), "contrast": lab,
                             "lfc": float(d.loc[g, "logfoldchanges"])})
    if rows:
        f = pd.DataFrame(rows)
        piv = f.pivot(index="gene", columns="contrast", values="lfc")
        piv = piv.reindex([g.replace("lncRNA:", "") for g in genes]).dropna(how="all")
        y = np.arange(len(piv))
        h = 0.36
        for k, (col, c) in enumerate(zip(piv.columns, [BLUE, RED])):
            axC.barh(y + (k - 0.5) * h, piv[col].values, height=h * 0.9,
                     color=c, label=col)
        axC.axvline(0, color=DARK, lw=0.8)
        axC.set_yticks(y)
        axC.set_yticklabels(piv.index, fontsize=8)
        axC.invert_yaxis()
        axC.set_xlabel("log$_2$FC in the nominal\ncocaine vs sucrose contrast",
                       fontsize=8.5)
        axC.legend(fontsize=7, frameon=False, loc="lower left")
    else:
        axC.text(0.5, 0.5, "run scripts/05 first", ha="center", va="center",
                 transform=axC.transAxes, fontsize=9, color=GREY)
    axC.set_title("C  Sex markers appear in\n     both treatment contrasts",
                  fontsize=10.5, loc="left")
    axC.spines[["top", "right"]].set_visible(False)
    axC.tick_params(labelsize=8)

    fig.suptitle("Sex-specific markers do not agree with the deposited sample labels",
                 fontsize=12.5, y=1.03)
    out = FIG_DIR / "fig_sexmarker_discrepancy.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out.name}")

    print("\nPer-sample means:")
    print(m.round(3).to_string())
    print("\nGrouped by declared treatment:")
    g = m.copy(); g["treatment"] = trt
    print(g.groupby("treatment").agg(["min", "max"]).round(3).to_string())


if __name__ == "__main__":
    main()
