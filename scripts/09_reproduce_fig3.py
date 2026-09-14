#!/usr/bin/env python3
"""
09_reproduce_fig3.py — reproduce Figure 3A-C of Baker et al. (2021).

WHAT THE ORIGINAL FIGURE SHOWS
  3A  Gene x cluster matrix for MALES: rows are clusters, columns are genes
      differentially expressed at |log_e FC| > 1.0 with Bonferroni-adjusted
      p < 0.05. Magenta = up in cocaine, turquoise = down. 133 genes.
  3B  The same for FEMALES. 54 genes.
  3C  Venn diagrams for six clusters the authors identify as having
      sexually divergent responses (their C10, C12, C15, C16, C19, C22),
      using the looser |log_e FC| > 0.5 threshold.

THRESHOLD NOTE
The paper uses natural-log fold changes; scanpy reports log2. Their
|log_e FC| > 1.0 is |log2FC| > 1.443, and |log_e FC| > 0.5 is |log2FC| > 0.721.
This script uses the converted values so the comparison is like-for-like,
rather than the |log2FC| > 1.0 used elsewhere in this analysis.

CLUSTER CORRESPONDENCE
The paper's cluster numbers are not ours. Panel C uses the cell types matched
to their divergent clusters by marker gene, where a match exists:

    their C16  -> Unannotated (antennal/optic lobe mix)
    their C22  -> Surface glia and fat body
    their C12  -> GABAergic neurons type 3
    their C10  -> Central brain B cholinergic neurons

Their C15 and C19 have no confident counterpart here. Rather than force a
match, panel C shows the four matched cell types plus the two cell types with
the largest female-male divergence in this dataset, each labelled with which
is which.

Run:  python scripts/09_reproduce_fig3.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Circle, Patch

from config import TABLE_DIR, FIG_DIR, PADJ_THRESHOLD

# Paper thresholds, converted from natural log to log2
LOG2_STRICT = 1.0 / np.log(2)   # their |log_e FC| > 1.0  -> 1.443
LOG2_LOOSE = 0.5 / np.log(2)    # their |log_e FC| > 0.5  -> 0.721

# Their divergent clusters mapped to our cell types where a match is defensible
PAPER_DIVERGENT = {
    "C10": "Central brain B cholinergic neurons",
    "C12": "GABAergic neurons type 3",
    "C16": "Unannotated (antennal/optic lobe mix)",
    "C22": "Surface glia and fat body",
}


def load_per_cluster():
    p = TABLE_DIR / "de_per_cluster_significant_genes.csv"
    if not p.exists():
        sys.exit(f"{p} not found — run scripts/05_de_analysis.py first.")
    df = pd.read_csv(p)
    df = df[~df["cluster"].astype(str).str.contains("excluded")]
    return df


def gene_cluster_matrix(df, sex, lfc, ax, title):
    """Panel A/B: binary gene x cluster matrix, coloured by direction."""
    d = df[(df.sex == sex) &
           (df.logfoldchanges.abs() > lfc) &
           (df.pvals_adj < PADJ_THRESHOLD)]
    if d.empty:
        ax.text(0.5, 0.5, f"No genes pass |log2FC| > {lfc:.2f}\nin {sex.lower()}s",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(title, fontsize=10)
        return 0, 0

    # +1 up in cocaine, -1 down, 0 not differentially expressed in that cluster
    mat = (d.assign(v=np.sign(d.logfoldchanges))
             .pivot_table(index="cluster", columns="names", values="v",
                          aggfunc="first")
             .fillna(0))

    # Order genes by how many clusters they appear in, then by direction, so
    # the block structure of the original figure is visible.
    order = (mat != 0).sum(axis=0).sort_values(ascending=False).index
    mat = mat[order]
    mat = mat.loc[(mat != 0).sum(axis=1).sort_values(ascending=False).index]

    cmap = ListedColormap(["#40C4C4", "#FFFFFF", "#D6219B"])  # down, none, up
    ax.imshow(mat.values, aspect="auto", cmap=cmap, vmin=-1, vmax=1,
              interpolation="nearest")
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels([c[:38] for c in mat.index], fontsize=6)
    ax.set_xticks([])
    ax.set_xlabel(f"{mat.shape[1]} differentially expressed genes", fontsize=8)
    ax.set_title(title, fontsize=10)
    for sp in ax.spines.values():
        sp.set_linewidth(0.4)
    return mat.shape[1], mat.shape[0]


def venn(ax, only_a, shared, only_b, label_a, label_b, title):
    """Two-circle Venn with counts. Areas are fixed, not proportional."""
    ax.set_xlim(-2.2, 2.2); ax.set_ylim(-1.6, 1.8)
    ax.set_aspect("equal"); ax.axis("off")

    ax.add_patch(Circle((-0.55, 0), 1.15, facecolor="#4FB39B", alpha=0.55,
                        edgecolor="black", lw=0.6))
    ax.add_patch(Circle((0.55, 0), 1.15, facecolor="#F2E85C", alpha=0.55,
                        edgecolor="black", lw=0.6))

    ax.text(-1.25, 0, str(only_a), ha="center", va="center", fontsize=10)
    ax.text(0, 0, str(shared), ha="center", va="center", fontsize=10)
    ax.text(1.25, 0, str(only_b), ha="center", va="center", fontsize=10)
    ax.text(-0.55, 1.35, label_a, ha="center", fontsize=8)
    ax.text(0.55, 1.35, label_b, ha="center", fontsize=8)
    ax.set_title(title, fontsize=8.5, pad=2)


def main():
    df = load_per_cluster()
    print(f"Loaded {len(df):,} significant per-cluster DE calls\n")

    # ---- Panels A and B --------------------------------------------------
    print("Panels A/B — gene x cluster matrices at the paper's strict threshold")
    print(f"  |log2FC| > {LOG2_STRICT:.3f}  (their |log_e FC| > 1.0)\n")

    fig, axes = plt.subplots(2, 1, figsize=(11, 9))
    ng_m, nc_m = gene_cluster_matrix(df, "Male", LOG2_STRICT, axes[0],
                                     "A  Male: cocaine-responsive genes by cell type")
    ng_f, nc_f = gene_cluster_matrix(df, "Female", LOG2_STRICT, axes[1],
                                     "B  Female: cocaine-responsive genes by cell type")
    print(f"  males:   {ng_m} genes across {nc_m} cell types "
          f"(Baker et al.: 133 genes)")
    print(f"  females: {ng_f} genes across {nc_f} cell types "
          f"(Baker et al.: 54 genes)")

    handles = [Patch(facecolor="#D6219B", edgecolor="black", lw=.4,
                     label="up in cocaine"),
               Patch(facecolor="#40C4C4", edgecolor="black", lw=.4,
                     label="down in cocaine")]
    axes[0].legend(handles=handles, loc="upper right", fontsize=7,
                   frameon=False, bbox_to_anchor=(1.0, 1.28), ncol=2)
    fig.suptitle("Reproduction of Baker et al. Figure 3A,B", fontsize=11, y=0.995)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig6ab_gene_cluster_matrix.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig6ab_gene_cluster_matrix.png")

    # ---- Panel C ---------------------------------------------------------
    print("\nPanel C — male/female overlap per cell type")
    print(f"  |log2FC| > {LOG2_LOOSE:.3f}  (their |log_e FC| > 0.5)\n")

    d = df[(df.logfoldchanges.abs() > LOG2_LOOSE) &
           (df.pvals_adj < PADJ_THRESHOLD)]

    sets = {}
    for ct in d.cluster.unique():
        m = set(d[(d.cluster == ct) & (d.sex == "Male")]["names"])
        f = set(d[(d.cluster == ct) & (d.sex == "Female")]["names"])
        if len(m | f) >= 5:
            sets[ct] = (f, m)

    # Matched cell types first, then whichever remaining types diverge most
    matched = [(p, ct) for p, ct in PAPER_DIVERGENT.items() if ct in sets]
    rest = sorted(
        (ct for ct in sets if ct not in dict(matched).values()),
        key=lambda c: -(len(sets[c][0] ^ sets[c][1])),
    )
    chosen = matched + [(None, ct) for ct in rest[:max(0, 6 - len(matched))]]

    rows = []
    n = len(chosen)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.1 * nrow))
    axes = np.atleast_1d(axes).ravel()

    for i, (paper_cl, ct) in enumerate(chosen):
        f, m = sets[ct]
        only_f, shared, only_m = len(f - m), len(f & m), len(m - f)
        tag = f"{ct[:30]}\n(≈ their {paper_cl})" if paper_cl else f"{ct[:30]}\n(no paper match)"
        venn(axes[i], only_f, shared, only_m, "♀", "♂", tag)
        rows.append({"cell_type": ct, "paper_cluster": paper_cl or "",
                     "female_only": only_f, "shared": shared, "male_only": only_m,
                     "jaccard": round(shared / max(1, len(f | m)), 3)})
        print(f"  {ct[:38]:<38} ♀only {only_f:>3}  shared {shared:>3}  ♂only {only_m:>3}")

    for j in range(n, len(axes)):
        axes[j].axis("off")

    fig.suptitle("Reproduction of Baker et al. Figure 3C — "
                 "male/female overlap of cocaine-responsive genes",
                 fontsize=11, y=1.0)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig6c_sex_overlap_venn.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig6c_sex_overlap_venn.png")

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "fig3c_sex_overlap.csv", index=False)
    print(f"\n{out.to_string(index=False)}")
    print("\nJaccard index is shared / union. The paper reports very little")
    print("overlap in its divergent clusters; compare directly.")


if __name__ == "__main__":
    main()
