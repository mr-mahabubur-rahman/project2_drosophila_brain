#!/usr/bin/env python3
"""
make_graphical_abstract.py — graphical abstract and analysis flowchart.

REWRITTEN after the sample labels were reverted to the deposited ones.

The earlier version of this script encoded a conclusion that no longer holds:
that the folder labels were scrambled, that the female response fell below its
own background, and that a male-specific mitochondrial confound drove the
result. Supplemental Table S2 showed the labels were correct as deposited, and
re-analysis under them gives different numbers and a different finding.

The current finding is that sex-specific markers do not agree with the deposited
sex labels, so neither treatment contrast can be separated from a sex contrast.
Both figures now carry that.

Numbers here are hard-coded. If the analysis is re-run and they change, edit
this file rather than the images.

Run:  python tools/make_graphical_abstract.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np

RED = "#C44E52"
BLUE = "#4C72B0"
GREEN = "#4E9A6A"
AMBER = "#D9843B"
GREY = "#8C8C8C"
LIGHT = "#F2F2F2"
DARK = "#2B2B2B"


def box(ax, x, y, w, h, text, fc="white", ec=DARK, lw=1.0, fs=8.5,
        bold=False, tc=DARK, align="center", pad=0.004, r=0.010):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad={pad},rounding_size={r}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2))
    if text:
        ax.text(x + w / 2 if align == "center" else x + 0.012,
                y + h / 2, text, ha=align, va="center", fontsize=fs,
                fontweight="bold" if bold else "normal",
                color=tc, zorder=3, linespacing=1.45)


def arrow(ax, x1, y1, x2, y2, color=GREY, lw=1.3, style="-|>", ls="-"):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle=style, mutation_scale=12,
        linewidth=lw, color=color, linestyle=ls,
        shrinkA=2, shrinkB=2, zorder=1))


# ======================================================================
# 1. GRAPHICAL ABSTRACT
# ======================================================================
def graphical_abstract(path):
    fig, ax = plt.subplots(figsize=(11.5, 5.9))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.970,
            "Independent reanalysis of the cocaine-exposed $\\it{Drosophila}$ brain",
            ha="center", fontsize=15, fontweight="bold", color=DARK)
    ax.text(0.5, 0.938,
            "The atlas reproduces. The sample labels do not agree with the data.",
            ha="center", fontsize=10.5, color=GREY, style="italic")

    TOP, BH = 0.688, 0.222

    # ---------------- input ----------------
    box(ax, 0.020, TOP, 0.205, BH, "", fc=LIGHT, ec=GREY, lw=0.8)
    ax.text(0.1225, TOP + 0.190, "DATA", ha="center", fontsize=8.5,
            fontweight="bold", color=GREY)
    ax.text(0.1225, TOP + 0.154, "GEO: GSE152495", ha="center", fontsize=10,
            fontweight="bold", color=DARK)
    ax.text(0.1225, TOP + 0.120, "8 samples, 88,991 cells", ha="center",
            fontsize=8.4, color=DARK)
    for i, (lbl, col) in enumerate([("\u2640 Sucrose", BLUE), ("\u2640 Cocaine", RED),
                                    ("\u2642 Sucrose", BLUE), ("\u2642 Cocaine", RED)]):
        x = 0.035 + (i % 2) * 0.094
        y = TOP + 0.068 - (i // 2) * 0.032
        ax.add_patch(Rectangle((x, y), 0.086, 0.025, facecolor=col,
                               alpha=0.22, edgecolor=col, linewidth=0.7))
        ax.text(x + 0.043, y + 0.0125, f"{lbl} \u00d72", ha="center", va="center",
                fontsize=7.2, color=DARK)

    arrow(ax, 0.233, TOP + BH / 2, 0.278, TOP + BH / 2, color=GREY, lw=1.8)

    # ---------------- pipeline ----------------
    box(ax, 0.288, TOP, 0.215, BH, "", fc=LIGHT, ec=GREY, lw=0.8)
    ax.text(0.3955, TOP + 0.190, "PIPELINE", ha="center", fontsize=8.5,
            fontweight="bold", color=GREY)
    ax.text(0.3955, TOP + 0.154, "Scanpy \u00b7 Leiden \u00b7 Wilcoxon", ha="center",
            fontsize=9.6, fontweight="bold", color=DARK)
    for i, line in enumerate([
            "3 data errors corrected",
            "86,177 cells retained",
            "(paper 86,224 \u2014 0.05% apart)",
            "30 clusters \u00b7 24 annotated"]):
        ax.text(0.3955, TOP + 0.116 - i * 0.029, line, ha="center",
                fontsize=8.2, color=DARK)

    arrow(ax, 0.511, TOP + BH / 2, 0.556, TOP + BH / 2, color=GREY, lw=1.8)

    # ---------------- outcomes ----------------
    box(ax, 0.566, TOP + 0.122, 0.414, 0.100, "", fc="#EAF3EE", ec=GREEN, lw=1.2)
    ax.text(0.583, TOP + 0.196, "THE ATLAS REPRODUCES", fontsize=8.6,
            fontweight="bold", color=GREEN, va="center")
    for i, line in enumerate(["cell count to 0.05% \u00b7 30 clusters \u00b7 major cell types",
                              "no batch effect \u00b7 Kenyon cells, surface glia recovered"]):
        ax.text(0.583, TOP + 0.166 - i * 0.027, line, fontsize=8.2,
                color=DARK, va="center")

    box(ax, 0.566, TOP, 0.414, 0.100, "", fc="#FBF0E6", ec=AMBER, lw=1.2)
    ax.text(0.583, TOP + 0.074, "THE LABELS DO NOT", fontsize=8.6,
            fontweight="bold", color=AMBER, va="center")
    for i, line in enumerate(["sex markers split the samples by TREATMENT,",
                              "crossing both declared sexes"]):
        ax.text(0.583, TOP + 0.044 - i * 0.027, line, fontsize=8.2,
                color=DARK, va="center")

    # ---------------- key result ----------------
    ax.plot([0.020, 0.980], [0.638, 0.638], color="#DDDDDD", lw=1)
    ax.text(0.5, 0.602, "KEY RESULT \u2014 sex markers do not match the sample labels",
            ha="center", fontsize=10.5, fontweight="bold", color=DARK)
    ax.text(0.5, 0.572,
            "Ranges across the four samples in each group. No sample bridges the gap.",
            ha="center", fontsize=8.4, color=GREY, style="italic")

    rows = [
        ("roX1", "male-specific", "3.85 \u2013 3.97", "0.06 \u2013 0.12"),
        ("roX2", "male-specific", "2.16 \u2013 2.40", "0.01 \u2013 0.02"),
        ("Yp1", "female-specific", "0.000 \u2013 0.001", "0.07 \u2013 0.26"),
        ("Yp3", "female-specific", "0.006 \u2013 0.010", "0.07 \u2013 0.43"),
    ]
    ytab = 0.520
    rh = 0.038
    x0, wg, wn, wv = 0.175, 0.115, 0.135, 0.180

    ax.text(x0 + wg / 2, ytab + 0.020, "marker", ha="center", fontsize=8.4,
            fontweight="bold", color=GREY)
    ax.text(x0 + wg + wn / 2, ytab + 0.020, "specificity", ha="center",
            fontsize=8.4, fontweight="bold", color=GREY)
    ax.text(x0 + wg + wn + wv / 2, ytab + 0.020, "4 Sucrose samples",
            ha="center", fontsize=8.6, fontweight="bold", color=BLUE)
    ax.text(x0 + wg + wn + wv * 1.5, ytab + 0.020, "4 Cocaine samples",
            ha="center", fontsize=8.6, fontweight="bold", color=RED)

    for i, (g, spec, suc, coc) in enumerate(rows):
        y = ytab - i * rh
        if i % 2 == 0:
            ax.add_patch(Rectangle((x0 - 0.01, y - rh * 0.42), wg + wn + wv * 2 + 0.02,
                                   rh * 0.84, facecolor="#F7F7F7", edgecolor="none"))
        ax.text(x0 + wg / 2, y, g, ha="center", va="center", fontsize=9,
                fontweight="bold", color=DARK, style="italic")
        ax.text(x0 + wg + wn / 2, y, spec, ha="center", va="center",
                fontsize=7.8, color=GREY)
        ax.text(x0 + wg + wn + wv / 2, y, suc, ha="center", va="center",
                fontsize=9, color=BLUE, fontweight="bold")
        ax.text(x0 + wg + wn + wv * 1.5, y, coc, ha="center", va="center",
                fontsize=9, color=RED, fontweight="bold")

    ax.text(0.5, ytab - 4 * rh - 0.012,
            "Yolk proteins are not transcribed in males. The highest value in the "
            "dataset is in a sample labelled male.",
            ha="center", fontsize=8.2, color=DARK, style="italic")

    # ---------------- consequence ----------------
    box(ax, 0.020, 0.212, 0.960, 0.098, "", fc="#FBF0E6", ec=AMBER, lw=1.0)
    ax.text(0.5, 0.284, "CONSEQUENCE", ha="center", fontsize=8.5,
            fontweight="bold", color=AMBER)
    ax.text(0.5, 0.257,
            "Comparing cocaine with sucrose within each declared sex returns the sex "
            "markers among the most significant genes in BOTH contrasts",
            ha="center", fontsize=8.6, color=DARK)
    ax.text(0.5, 0.230,
            "roX1 and roX2 at approximately \u22129.5 log$_2$FC, adjusted $p$ below "
            "machine precision",
            ha="center", fontsize=8.6, color=DARK)

    # ---------------- verification + conclusion ----------------
    box(ax, 0.020, 0.088, 0.960, 0.108, "", fc="#F7F7F7", ec=GREY, lw=0.9)
    ax.text(0.5, 0.172,
            "Sample identity verified three ways, all in agreement: "
            "Supplemental Table S2 cell counts (8/8), GEO titles, authors' code",
            ha="center", fontsize=8.4, color=DARK)
    ax.text(0.5, 0.146,
            "The discrepancy did not arise during analysis or deposition.",
            ha="center", fontsize=8.4, color=DARK)
    ax.text(0.5, 0.112,
            "No differential expression result from these data can presently be "
            "attributed to cocaine. Referred to the original authors.",
            ha="center", fontsize=9.0, fontweight="bold", color=DARK)

    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {path}")


# ======================================================================
# 2. FLOWCHART
# ======================================================================
def flowchart(path):
    fig, ax = plt.subplots(figsize=(8.6, 14.2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.986, "Analysis workflow", ha="center", fontsize=15,
            fontweight="bold", color=DARK)
    ax.text(0.5, 0.968, "Bracketed numbers are the pipeline scripts",
            ha="center", fontsize=8.2, color=GREY, style="italic")

    CX, W, GAP = 0.50, 0.66, 0.014
    y = 0.955

    steps = [
        ("GSE152495 \u2014 8 CellRanger matrices",
         ["88,991 cells \u00d7 17,481 genes"], LIGHT, GREY),
        ("Pre-analysis corrections  [tools/]",
         ["space-delimited features.tsv \u00b7 gene named 'nan'",
          "five vertebrate ortholog symbols (trh = trachealess)"],
         "#FBF0E6", AMBER),
        ("Sample identity verified  [tools/]",
         ["Table S2 cell counts 8/8 \u00b7 GEO titles \u00b7 authors' code",
          "labels used AS DEPOSITED"], "#EAF3EE", GREEN),
        ("Load and merge  [01]",
         ["metadata attached to every cell"], "white", GREY),
        ("Quality control  [02]",
         ["300\u20132500 genes \u00b7 <10% mito \u00b7 \u22655 cells per gene",
          "86,177 cells retained (paper: 86,224)"], "white", GREY),
        ("Normalise, HVG, PCA, cluster  [03]",
         ["CPM+log1p \u00b7 2,000 HVG \u00b7 30 PCs \u00b7 Leiden 0.8",
          "30 clusters (paper: 36); sweep 0.4\u20131.4"], "white", GREY),
        ("Cell-type annotation  [04]",
         ["two independent marker panels",
          "24 of 30 assigned; 1 cluster excluded"], "white", GREY),
        ("Sex-marker check  [tools/]",
         ["roX1/roX2 and Yp1\u2013Yp3 vs declared labels",
          "MARKERS SPLIT BY TREATMENT, NOT SEX"], "#FBF0E6", AMBER),
    ]

    centres = []
    for title, lines, fc, ec in steps:
        h = 0.030 + 0.017 * len(lines)
        box(ax, CX - W / 2, y - h, W, h, "", fc=fc, ec=ec, lw=1.0)
        ax.text(CX, y - 0.016, title, ha="center", fontsize=9.4,
                fontweight="bold", color=DARK, va="center")
        for i, ln in enumerate(lines):
            ax.text(CX, y - 0.034 - i * 0.017, ln, ha="center", fontsize=7.7,
                    color=DARK, va="center")
        centres.append((y, y - h))
        y = y - h - GAP

    for i in range(len(centres) - 1):
        arrow(ax, CX, centres[i][1], CX, centres[i + 1][0], color=GREY, lw=1.5)

    # ---- controls ----
    y_ctrl_top = y - 0.006
    ax.plot([0.02, 0.98], [y_ctrl_top, y_ctrl_top], color="#DDDDDD", lw=1)
    ax.text(0.5, y_ctrl_top - 0.026, "DIFFERENTIAL EXPRESSION AND CONTROLS",
            ha="center", fontsize=9.6, fontweight="bold", color=BLUE)
    ax.text(0.5, y_ctrl_top - 0.045,
            "reported for completeness; each contrast is confounded with sex",
            ha="center", fontsize=7.9, color=GREY, style="italic")
    arrow(ax, CX, centres[-1][1], CX, y_ctrl_top, color=GREY, lw=1.5)

    controls = [
        ("Differential expression  [05]", "Cocaine vs sucrose within\neach declared sex",
         "90 and 152 genes \u2014 but\nsex markers top both lists", False, "Fig 3C"),
        ("Pseudobulk  [05b]", "Do cell-level results hold\nat replicate level?",
         "100% direction agreement\n\u03c1 = 0.94, 0.96", True, "Fig 17"),
        ("Interaction model  [07]", "Does any gene respond\ndifferently by sex?",
         "0 genes at FDR < 0.05\n(4 residual df)", False, "Fig 14"),
        ("Cluster-size matching  [08B]", "Is the ranking power\nor biology?",
         "\u03c1 = 0.63 size vs DE count;\nranking shifts when matched", False, "Fig 15"),
        ("Permutation control  [10]", "Is the effect larger than\nits own background?",
         "Male 90 vs 16 / 44\nFemale 152 vs 19 / 33", True, "Fig 16"),
    ]

    cy = y_ctrl_top - 0.054
    CH, CGAP = 0.038, 0.008
    for name, question, result, good, figref in controls:
        top, bot = cy, cy - CH
        box(ax, 0.045, bot, 0.345, CH, "", fc="#EEF2F8", ec=BLUE, lw=0.9)
        ax.text(0.058, top - 0.012, name, fontsize=8.3, fontweight="bold",
                color=BLUE, va="center")
        for i, ln in enumerate(question.split("\n")):
            ax.text(0.058, top - 0.025 - i * 0.0125, ln, fontsize=7.0,
                    color=DARK, va="center")
        fcol = "#EAF3EE" if good else "#FBF0E6"
        ecol = GREEN if good else AMBER
        box(ax, 0.455, bot, 0.345, CH, "", fc=fcol, ec=ecol, lw=0.9)
        for i, ln in enumerate(result.split("\n")):
            ax.text(0.468, top - 0.017 - i * 0.016, ln, fontsize=7.4,
                    color=DARK, va="center")
        arrow(ax, 0.393, top - CH / 2, 0.452, top - CH / 2, color=GREY, lw=1.1)
        ax.text(0.815, top - CH / 2, figref, fontsize=7.4, color=GREY,
                va="center", style="italic")
        cy = bot - CGAP

    # ---- conclusion ----
    ctop = cy - 0.006
    CHH = 0.070
    box(ax, 0.045, ctop - CHH, 0.755, CHH, "", fc="#F7F7F7", ec=DARK, lw=1.1)
    ax.text(0.4225, ctop - 0.014, "CONCLUSION", ha="center", fontsize=8.4,
            fontweight="bold", color=DARK)
    ax.text(0.4225, ctop - 0.032,
            "The atlas reproduces: 86,177 cells, 30 clusters, major cell types recovered.",
            ha="center", fontsize=8.0, color=DARK)
    ax.text(0.4225, ctop - 0.047,
            "Sex markers do not agree with the deposited labels, so neither treatment",
            ha="center", fontsize=8.0, color=DARK)
    ax.text(0.4225, ctop - 0.062,
            "contrast can be separated from a sex contrast. Referred to the authors.",
            ha="center", fontsize=8.0, color=DARK)

    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    graphical_abstract("results/figures/fig0_graphical_abstract.png")
    flowchart("results/figures/fig0b_flowchart.png")
