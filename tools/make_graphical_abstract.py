#!/usr/bin/env python3
"""
make_graphical_abstract.py — graphical abstract and analysis flowchart.

Both are drawn rather than screenshotted so they can be regenerated if numbers
change. Colour is used consistently across the two figures and the rest of the
report: red for the treatment/real contrast, blue for nulls and controls, grey
for neutral material, green for findings that reproduce and amber for those
that do not survive testing.
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
    ax.text(x + w / 2 if align == "center" else x + 0.012,
            y + h / 2, text,
            ha=align, va="center", fontsize=fs,
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
    fig, ax = plt.subplots(figsize=(11.5, 6.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.968,
            "Independent reanalysis of the cocaine-exposed $\\it{Drosophila}$ brain",
            ha="center", fontsize=15, fontweight="bold", color=DARK)
    ax.text(0.5, 0.932,
            "Which published findings survive explicit statistical controls?",
            ha="center", fontsize=10.5, color=GREY, style="italic")

    TOP, BH = 0.635, 0.245          # box bottom and height

    # ---------------- input ----------------
    box(ax, 0.020, TOP, 0.205, BH, "", fc=LIGHT, ec=GREY, lw=0.8)
    ax.text(0.1225, TOP + 0.212, "DATA", ha="center", fontsize=8.5,
            fontweight="bold", color=GREY)
    ax.text(0.1225, TOP + 0.176, "GEO: GSE152495", ha="center", fontsize=10,
            fontweight="bold", color=DARK)
    ax.text(0.1225, TOP + 0.140, "8 samples, 88,991 cells", ha="center",
            fontsize=8.4, color=DARK)
    for i, (lbl, col) in enumerate([("\u2640 Sucrose", BLUE), ("\u2640 Cocaine", RED),
                                    ("\u2642 Sucrose", BLUE), ("\u2642 Cocaine", RED)]):
        x = 0.035 + (i % 2) * 0.094
        y = TOP + 0.085 - (i // 2) * 0.034
        ax.add_patch(Rectangle((x, y), 0.086, 0.026, facecolor=col,
                               alpha=0.22, edgecolor=col, linewidth=0.7))
        ax.text(x + 0.043, y + 0.013, f"{lbl} \u00d72", ha="center", va="center",
                fontsize=7.2, color=DARK)

    arrow(ax, 0.233, TOP + BH / 2, 0.278, TOP + BH / 2, color=GREY, lw=1.8)

    # ---------------- pipeline ----------------
    box(ax, 0.288, TOP, 0.215, BH, "", fc=LIGHT, ec=GREY, lw=0.8)
    ax.text(0.3955, TOP + 0.212, "PIPELINE", ha="center", fontsize=8.5,
            fontweight="bold", color=GREY)
    ax.text(0.3955, TOP + 0.176, "Scanpy \u00b7 Leiden \u00b7 Wilcoxon", ha="center",
            fontsize=9.6, fontweight="bold", color=DARK)
    for i, line in enumerate([
            "3 data errors corrected",
            "86,177 cells retained",
            "(paper 86,224 \u2014 0.05% apart)",
            "30 clusters \u00b7 24 annotated"]):
        ax.text(0.3955, TOP + 0.132 - i * 0.031, line, ha="center",
                fontsize=8.2, color=DARK)

    arrow(ax, 0.511, TOP + BH / 2, 0.556, TOP + BH / 2, color=GREY, lw=1.8)

    # ---------------- outcomes ----------------
    box(ax, 0.566, TOP + 0.136, 0.414, 0.109, "", fc="#EAF3EE", ec=GREEN, lw=1.2)
    ax.text(0.583, TOP + 0.218, "REPRODUCES", fontsize=8.6, fontweight="bold",
            color=GREEN, va="center")
    for i, line in enumerate(["cell count \u00b7 atlas structure \u00b7 no batch effect",
                              "near-zero overlap between sexes \u00b7 male response"]):
        ax.text(0.583, TOP + 0.186 - i * 0.028, line, fontsize=8.3,
                color=DARK, va="center")

    box(ax, 0.566, TOP, 0.414, 0.109, "", fc="#FBF0E6", ec=AMBER, lw=1.2)
    ax.text(0.583, TOP + 0.082, "DOES NOT SURVIVE TESTING", fontsize=8.6,
            fontweight="bold", color=AMBER, va="center")
    for i, line in enumerate(["female response \u2014 below its own background",
                              "male:female ratio \u2014 0.97\u00d7 to 6.43\u00d7 by metric"]):
        ax.text(0.583, TOP + 0.050 - i * 0.028, line, fontsize=8.3,
                color=DARK, va="center")

    # ---------------- key result panel ----------------
    ax.plot([0.020, 0.980], [0.585, 0.585], color="#DDDDDD", lw=1)
    ax.text(0.5, 0.548, "KEY RESULT \u2014 treatment effect against its own background",
            ha="center", fontsize=10.5, fontweight="bold", color=DARK)
    ax.text(0.5, 0.516,
            "Four samples per sex permit two null contrasts besides the true treatment split",
            ha="center", fontsize=8.4, color=GREY, style="italic")

    BASE, TOPBAR = 0.245, 0.462       # plotting area for bars
    for pi, (sex, vals, note) in enumerate([
            ("Male", [90, 7, 14], "treatment dominates"),
            ("Female", [14, 39, 35], "replicate axis dominates")]):
        x0 = 0.075 + pi * 0.48
        w = 0.36
        ax.text(x0 + w / 2, TOPBAR + 0.026, sex, ha="center", fontsize=11,
                fontweight="bold", color=DARK)
        ax.text(x0 + w / 2, TOPBAR + 0.002, note, ha="center", fontsize=8.2,
                color=GREY, style="italic")
        vmax, bw = 100.0, 0.075
        for bi, (v, c, lab) in enumerate(zip(
                vals, [RED, BLUE, GREY],
                ["treatment\n(real)", "replicate\n(null)", "diagonal\n(null)"])):
            bx = x0 + 0.045 + bi * 0.104
            bh = (v / vmax) * (TOPBAR - BASE - 0.038)
            ax.add_patch(Rectangle((bx, BASE), bw, bh, facecolor=c,
                                   edgecolor="none", zorder=2))
            ax.text(bx + bw / 2, BASE + bh + 0.013, str(v), ha="center",
                    fontsize=10, fontweight="bold", color=DARK)
            ax.text(bx + bw / 2, BASE - 0.020, lab, ha="center", fontsize=7.6,
                    color=DARK, linespacing=1.5, va="top")
        ax.plot([x0 + 0.032, x0 + w], [BASE, BASE], color=DARK, lw=1.1)
        ax.text(x0 + 0.020, BASE + (TOPBAR - BASE - 0.038) / 2, "DE genes",
                rotation=90, va="center", ha="center", fontsize=8.2, color=GREY)

    # ---------------- bottom conclusion ----------------
    box(ax, 0.020, 0.016, 0.960, 0.098, "", fc="#F7F7F7", ec=GREY, lw=0.9)
    ax.text(0.5, 0.086,
            "In males the treatment effect exceeds both nulls by 6\u201313\u00d7.",
            ha="center", fontsize=8.8, color=DARK)
    ax.text(0.5, 0.060,
            "In females both nulls exceed it, and 64% of the 14 female genes also "
            "appear in a null contrast.",
            ha="center", fontsize=8.8, color=DARK)
    ax.text(0.5, 0.031,
            "Design-level conclusions reproduce; the female response is not "
            "resolvable above between-sample variation.",
            ha="center", fontsize=9.2, fontweight="bold", color=DARK)

    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {path}")


# ======================================================================
# 2. FLOWCHART
# ======================================================================
def flowchart(path):
    fig, ax = plt.subplots(figsize=(8.6, 13.0))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.986, "Analysis workflow", ha="center", fontsize=15,
            fontweight="bold", color=DARK)
    ax.text(0.5, 0.968,
            "Bracketed numbers are the pipeline scripts",
            ha="center", fontsize=8.2, color=GREY, style="italic")

    CX, W = 0.50, 0.66
    GAP = 0.016
    y = 0.952                       # running TOP edge

    steps = [
        ("GSE152495 \u2014 8 CellRanger matrices",
         ["88,991 cells \u00d7 17,481 genes"], LIGHT, GREY),
        ("Pre-analysis corrections  [tools/]",
         ["delimiter \u00b7 gene named 'nan' \u00b7 scrambled sample labels",
          "sex confirmed by markers; treatment inferred from GEO order"],
         "#FBF0E6", AMBER),
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
        ("Differential expression  [05]",
         ["Wilcoxon, cocaine vs sucrose",
          "pooled \u00b7 per cluster \u00b7 sex overlap"], "white", GREY),
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
    ax.text(0.5, y_ctrl_top - 0.026, "STATISTICAL CONTROLS", ha="center",
            fontsize=10, fontweight="bold", color=BLUE)
    ax.text(0.5, y_ctrl_top - 0.045,
            "each answers a question the gene count alone cannot",
            ha="center", fontsize=7.9, color=GREY, style="italic")
    arrow(ax, CX, centres[-1][1], CX, y_ctrl_top, color=GREY, lw=1.5)

    controls = [
        ("Pseudobulk  [05b]", "Do cell-level results hold\nat replicate level?",
         "100% direction agreement\n\u03c1 = 0.78 (male)", True, "Fig 17"),
        ("Interaction model  [07]", "Does any gene respond\ndifferently by sex?",
         "0 genes at FDR < 0.05\n(4 residual df)", False, "Fig 14"),
        ("Gene-class exclusion  [08A]", "Do mitochondrial genes\ndrive the male bias?",
         "No \u2014 ratio unchanged\n(6.43\u00d7 \u2192 6.83\u00d7)", False, "Fig 11"),
        ("Cluster-size matching  [08B]", "Is the ranking power\nor biology?",
         "Kenyon 1st \u2192 7th;\nsurface glia stable at 3rd", False, "Fig 15"),
        ("Permutation control  [10]", "Is the effect larger than\nits own background?",
         "Male 90 vs 7 / 14  \u2713\nFemale 14 vs 39 / 35  \u2717", False, "Fig 16"),
    ]

    cy = y_ctrl_top - 0.058
    CH = 0.046
    CGAP = 0.011
    for name, question, result, good, figref in controls:
        top = cy
        bot = cy - CH
        box(ax, 0.045, bot, 0.345, CH, "", fc="#EEF2F8", ec=BLUE, lw=0.9)
        ax.text(0.058, top - 0.012, name, fontsize=8.3, fontweight="bold",
                color=BLUE, va="center")
        for i, ln in enumerate(question.split("\\n")):
            ax.text(0.058, top - 0.027 - i * 0.014, ln, fontsize=7.2,
                    color=DARK, va="center")
        fcol = "#EAF3EE" if good else "#FBF0E6"
        ecol = GREEN if good else AMBER
        box(ax, 0.455, bot, 0.345, CH, "", fc=fcol, ec=ecol, lw=0.9)
        for i, ln in enumerate(result.split("\\n")):
            ax.text(0.468, top - 0.017 - i * 0.016, ln, fontsize=7.4,
                    color=DARK, va="center")
        arrow(ax, 0.393, top - CH / 2, 0.452, top - CH / 2, color=GREY, lw=1.1)
        ax.text(0.815, top - CH / 2, figref, fontsize=7.4, color=GREY,
                va="center", style="italic")
        cy = bot - CGAP

    # ---- conclusion ----
    ctop = cy - 0.010
    CHH = 0.066
    box(ax, 0.045, ctop - CHH, 0.755, CHH, "", fc="#F7F7F7", ec=DARK, lw=1.1)
    ax.text(0.4225, ctop - 0.016, "CONCLUSION", ha="center", fontsize=8.4,
            fontweight="bold", color=DARK)
    ax.text(0.4225, ctop - 0.035,
            "Design-level findings reproduce. Gene-level findings are analysis-dependent,",
            ha="center", fontsize=8.1, color=DARK)
    ax.text(0.4225, ctop - 0.052,
            "and the female response is not resolvable above between-sample variation.",
            ha="center", fontsize=8.1, color=DARK)
    arrow(ax, 0.4225, cy + CGAP - CH - 0.001, 0.4225, ctop, color=GREY, lw=1.4)

    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    graphical_abstract("results/figures/fig0_graphical_abstract.png")
    flowchart("results/figures/fig0b_flowchart.png")
