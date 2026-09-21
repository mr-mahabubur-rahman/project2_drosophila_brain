#!/usr/bin/env python3
"""
make_sexmarker_figure.py — the labelling discrepancy in one figure.

WHAT IT SHOWS
Panel A: per-sample mean expression of the male-specific markers (roX1, roX2;
         top) and the female-specific yolk proteins (Yp1-Yp3; bottom), with
         samples ordered by declared sex. If the labels were correct, the male
         markers would be high in the four male-labelled samples. They are not:
         the split follows treatment.
Panel B: the same data as a two-marker scatter (roX1 against Yp1). Each point is
         one sample, coloured by declared treatment and shaped by declared sex.
         Samples separate by colour, not by shape.
Panel C: log2 fold changes of the five markers in the two nominal cocaine-vs-
         sucrose contrasts, showing that both contrasts return them.

LAYOUT CHANGES (v2)
- Gene identity in panel A is shown by shade, with a legend, instead of rotated
  text above every bar (which collided with the declared-sex labels).
- Declared-sex labels sit above panel A, clear of the bars.
- Panel B keeps log axes (so all eight samples are visible) but floors exact
  zeros at 1e-4 and says so on the plot. Group annotations sit clear of the
  points.
- The note under panel C is computed from the DE tables: it names the genes
  whose adjusted p is 0 in both contrasts and gives the largest remaining one.
- Per-sample means are also written to results/tables/sexmarker_per_sample_means.csv.

Panel C reads from the differential expression tables, so run scripts/05 first.

Run:  python tools/make_sexmarker_figure.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from config import H5AD_QC, FIG_DIR, TABLE_DIR

sc.settings.verbosity = 0

RED, BLUE, DARK, GREY = "#C44E52", "#4C72B0", "#2B2B2B", "#8C8C8C"

MALE_CANDIDATES = [("lncRNA:roX1", "roX1"), ("lncRNA:roX2", "roX2")]
FEMALE_CANDIDATES = ["Yp1", "Yp2", "Yp3"]
SUPPORT_GENES = ["Sxl", "msl-2", "mle", "tra"]   # tabulated, not plotted

ORDER = ["Female_Sucrose_R1", "Female_Sucrose_R2",
         "Female_Cocaine_R1", "Female_Cocaine_R2",
         "Male_Sucrose_R1", "Male_Sucrose_R2",
         "Male_Cocaine_R1", "Male_Cocaine_R2"]
SHORT = ["F Suc 1", "F Suc 2", "F Coc 1", "F Coc 2",
         "M Suc 1", "M Suc 2", "M Coc 1", "M Coc 2"]

MALE_ALPHA = [1.0, 0.50]
FEMALE_ALPHA = [1.0, 0.62, 0.34]


def resolve(adata):
    male = []
    for full, short in MALE_CANDIDATES:
        if full in adata.var_names:
            male.append(full)
        elif short in adata.var_names:
            male.append(short)
    female = [g for g in FEMALE_CANDIDATES if g in adata.var_names]
    if len(male) < 2 or len(female) < 3:
        found = [v for v in adata.var_names if "roX" in v or v.startswith("Yp")]
        sys.exit(f"markers not found; available roX/Yp: {found}")
    support = [g for g in SUPPORT_GENES if g in adata.var_names]
    return male, female, support


def clean(g):
    return g.replace("lncRNA:", "")


def bars(ax, m, genes, alphas, colors):
    x = np.arange(len(m))
    w = 0.8 / len(genes)
    for i, (g, a) in enumerate(zip(genes, alphas)):
        ax.bar(x + i * w - 0.4 + w / 2, m[g].values, width=w * 0.92,
               color=colors, alpha=a, edgecolor="none")
    ax.axvline(3.5, color=DARK, lw=1.0, ls="--")
    ax.set_xlim(-0.6, len(m) - 0.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    return x


def main():
    adata = sc.read_h5ad(H5AD_QC)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    male_g, female_g, support_g = resolve(adata)
    genes = male_g + female_g

    d = sc.get.obs_df(adata, keys=genes + support_g)
    d["sample"] = adata.obs["sample"].astype(str).values
    m = d.groupby("sample", observed=True).mean()
    order = [s for s in ORDER if s in m.index]
    if len(order) != len(m):
        sys.exit(f"sample names differ from ORDER: {sorted(m.index)}")
    short = [SHORT[ORDER.index(s)] for s in order]
    m = m.loc[order]
    trt = ["Cocaine" if "Cocaine" in s else "Sucrose" for s in order]
    colors = [RED if t == "Cocaine" else BLUE for t in trt]

    fig = plt.figure(figsize=(13.5, 5.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.45, 0.95, 1.05],
                          height_ratios=[1, 1], wspace=0.45, hspace=0.18)

    # ---------------- Panel A: grouped bars, two rows ----------------
    axA1 = fig.add_subplot(gs[0, 0])
    axA2 = fig.add_subplot(gs[1, 0], sharex=axA1)

    x = bars(axA1, m, male_g, MALE_ALPHA, colors)
    axA1.set_ylabel("Male markers\n(mean log-norm.)", fontsize=8.5)
    axA1.set_ylim(0, m[male_g].values.max() * 1.12)
    plt.setp(axA1.get_xticklabels(), visible=False)
    axA1.legend(handles=[Patch(facecolor=GREY, alpha=a, label=clean(g))
                         for g, a in zip(male_g, MALE_ALPHA)],
                fontsize=7, frameon=False, loc="upper right", ncol=2,
                handlelength=1.2, columnspacing=0.8)

    bars(axA2, m, female_g, FEMALE_ALPHA, colors)
    axA2.set_ylabel("Female markers\n(mean log-norm.)", fontsize=8.5)
    axA2.set_ylim(0, m[female_g].values.max() * 1.15)
    axA2.set_xticks(x)
    axA2.set_xticklabels(short, rotation=45, ha="right", fontsize=8)
    axA2.legend(handles=[Patch(facecolor=GREY, alpha=a, label=g)
                         for g, a in zip(female_g, FEMALE_ALPHA)],
                fontsize=7, frameon=False, loc="upper left", ncol=3,
                handlelength=1.2, columnspacing=0.8)

    tr = mtransforms.blended_transform_factory(axA1.transData, axA1.transAxes)
    axA1.text(1.5, 1.03, "declared female", transform=tr, ha="center",
              va="bottom", fontsize=8.5, color=DARK)
    axA1.text(5.5, 1.03, "declared male", transform=tr, ha="center",
              va="bottom", fontsize=8.5, color=DARK)
    axA1.set_title("A  Sex markers by sample", fontsize=10.5, loc="left", pad=20)

    fig.legend(handles=[Patch(facecolor=BLUE, label="declared sucrose"),
                        Patch(facecolor=RED, label="declared cocaine")],
               fontsize=7.5, frameon=False, ncol=2, loc="lower center",
               bbox_to_anchor=(0.25, -0.10))

    # ---------------- Panel B: linear two-marker scatter ----------------
    axB = fig.add_subplot(gs[:, 1])
    xg, yg = male_g[0], female_g[0]
    FLOOR = 1e-4                       # log axes cannot show an exact zero
    xv = m[xg].clip(lower=FLOOR)
    yv = m[yg].clip(lower=FLOOR)
    floored = bool((m[[xg, yg]] < FLOOR).any().any())
    for s, c in zip(order, colors):
        sex = "Female" if s.startswith("Female") else "Male"
        axB.scatter(xv[s], yv[s], s=110, c=c,
                    marker="o" if sex == "Female" else "^",
                    edgecolor=DARK, linewidth=0.7, zorder=3, alpha=0.9)
    axB.set_xscale("log"); axB.set_yscale("log")
    axB.set_xlim(xv.min() / 2.5, xv.max() * 2.5)
    axB.set_ylim(yv.min() / 4, yv.max() * 4)
    axB.grid(True, which="major", color="#EEEEEE", lw=0.6, zorder=0)

    coc = [s for s, t in zip(order, trt) if t == "Cocaine"]
    suc = [s for s, t in zip(order, trt) if t == "Sucrose"]
    gmean = lambda v: float(np.exp(np.log(v).mean()))
    axB.text(xv[coc].max() * 2.2, gmean(yv[coc]),
             "4 cocaine-labelled\n(2 declared F,\n 2 declared M)",
             fontsize=7.8, color=RED, ha="left", va="center")
    axB.text(xv[suc].min() / 2.2, gmean(yv[suc]),
             "4 sucrose-labelled\n(2 declared F,\n 2 declared M)",
             fontsize=7.8, color=BLUE, ha="right", va="center")
    if floored:
        axB.text(0.99, 0.01, f"values < {FLOOR:g} shown at {FLOOR:g}",
                 transform=axB.transAxes, ha="right", va="bottom",
                 fontsize=6.5, color=GREY, style="italic")

    axB.set_xlabel(f"{clean(xg)}  (male-specific)", fontsize=8.5)
    axB.set_ylabel(f"{yg}  (female-specific)", fontsize=8.5)
    axB.set_title("B  Samples separate by declared\n"
                  "    treatment, not declared sex", fontsize=10.5, loc="left")
    axB.spines[["top", "right"]].set_visible(False)
    axB.tick_params(labelsize=8)
    axB.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GREY,
               markeredgecolor=DARK, markersize=8, label="declared female"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=GREY,
               markeredgecolor=DARK, markersize=8, label="declared male"),
    ], fontsize=7.5, frameon=False, loc="upper center", ncol=2,
        bbox_to_anchor=(0.5, -0.12))

    # ---------------- Panel C: fold changes in both contrasts ----------------
    axC = fig.add_subplot(gs[:, 2])
    rows, padj = [], {}
    for sex, lab in [("female", "female-labelled"), ("male", "male-labelled")]:
        p = TABLE_DIR / f"de_global_{sex}_all.csv"
        if not p.exists():
            continue
        dd = pd.read_csv(p).set_index("names")
        for g in genes:
            if g in dd.index:
                rows.append({"gene": clean(g), "contrast": lab,
                             "lfc": float(dd.loc[g, "logfoldchanges"])})
                padj.setdefault(clean(g), []).append(float(dd.loc[g, "pvals_adj"]))
    if rows:
        f = pd.DataFrame(rows)
        piv = f.pivot(index="gene", columns="contrast", values="lfc")
        piv = piv.reindex([clean(g) for g in genes]).dropna(how="all")
        piv = piv[[c for c in ["female-labelled", "male-labelled"] if c in piv.columns]]
        y = np.arange(len(piv))
        h = 0.36
        for k, col in enumerate(piv.columns):
            c = BLUE if col == "female-labelled" else RED
            axC.barh(y + (k - 0.5) * h, piv[col].values, height=h * 0.9,
                     color=c, label=col)
        axC.axvline(0, color=DARK, lw=0.8)
        axC.set_yticks(y)
        axC.set_yticklabels(piv.index, fontsize=8.5, style="italic")
        axC.invert_yaxis()
        axC.set_xlabel("log$_2$FC in the nominal\ncocaine vs sucrose contrast",
                       fontsize=8.5)
        axC.legend(fontsize=7.5, frameon=False, loc="lower left")
        zero = [g for g, v in padj.items() if max(v) == 0.0]
        rest = [p for g, v in padj.items() if g not in zero for p in v]
        note = ("adj. p = 0 in both contrasts: "
                + ", ".join(zero)) if zero else ""
        if rest:
            note += ("\n" if note else "") + f"all others: adj. p \u2264 {max(rest):.2g}"
        axC.text(0.5, -0.26, note, transform=axC.transAxes, ha="center",
                 va="top", fontsize=7, color=GREY, style="italic")
    else:
        axC.text(0.5, 0.5, "run scripts/05 first", ha="center", va="center",
                 transform=axC.transAxes, fontsize=9, color=GREY)
    axC.set_title("C  Sex markers are significant\n"
                  "    in BOTH treatment contrasts",
                  fontsize=10.5, loc="left")
    axC.spines[["top", "right"]].set_visible(False)
    axC.tick_params(labelsize=8)

    fig.suptitle("Sex-specific markers do not agree with the deposited sample labels",
                 fontsize=12.5, y=1.04)
    out = FIG_DIR / "fig_sexmarker_discrepancy.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out.name}")

    tab = m.copy()
    tab.columns = [clean(c) for c in tab.columns]
    tab.insert(0, "declared_treatment", trt)
    tab.insert(0, "declared_sex", ["Female" if s.startswith("Female") else "Male" for s in order])
    tab.round(4).to_csv(TABLE_DIR / "sexmarker_per_sample_means.csv")
    print("wrote sexmarker_per_sample_means.csv")

    print("\nPer-sample means:")
    print(tab.round(3).to_string())
    print("\nRange by declared treatment:")
    print(tab.drop(columns="declared_sex").groupby("declared_treatment")
             .agg(["min", "max"]).round(3).T.to_string())
    print("\nadjusted p (female-labelled, male-labelled):")
    for g, v in padj.items():
        print(f"  {g:6s} " + "  ".join(f"{x:.3g}" for x in v))


if __name__ == "__main__":
    main()
