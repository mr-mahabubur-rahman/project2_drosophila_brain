#!/usr/bin/env python3
"""
rebuild_figures.py — regenerate every report figure from current checkpoints.

WHY THIS EXISTS
The sample labels were reverted to the deposited ones after Supplemental Table
S2 confirmed them. Every figure that depends on sex or treatment assignment was
therefore built on superseded data and must be redrawn. Clustering and
annotation figures are unaffected in content but are regenerated anyway so that
every file in results/figures/ has a consistent provenance.

Two figures (the DE-burden heatmap and the ratio comparison) were previously
produced by ad-hoc code typed at the prompt. They are included here so the
figure set is reproducible from scripts alone.

It also deletes the scanpy default plots that scripts 02 and 04 rewrite on every
run. Those are superseded by the fig1/fig3 versions and only cause confusion
when both are present.

Run:  python tools/rebuild_figures.py
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config

RED, BLUE, PURPLE, GREY = "#C44E52", "#4C72B0", "#8172B2", "#999999"


def run(script):
    """Run one of the standalone figure tools."""
    path = ROOT / "tools" / script
    if not path.exists():
        print(f"  SKIP {script} (not found)")
        return
    print(f"  running {script}")
    r = subprocess.run([sys.executable, str(path)], capture_output=True, text=True)
    if r.returncode != 0:
        print(f"    FAILED:\n{r.stdout[-800:]}\n{r.stderr[-800:]}")
    else:
        for line in r.stdout.strip().splitlines()[-4:]:
            print(f"    {line}")


def de_burden_heatmap():
    """Figure: DE genes per cell type and sex, excluded cluster omitted."""
    s = pd.read_csv(config.TABLE_DIR / "de_per_cluster_summary.csv")
    s = s.dropna(subset=["n_sig"])
    s = s[~s.cluster.astype(str).str.contains("excluded")]
    p = s.pivot_table(index="cluster", columns="sex", values="n_sig").fillna(0)
    p = p.loc[p.sum(axis=1).sort_values().index]

    fig, ax = plt.subplots(figsize=(7, 8))
    im = ax.imshow(p.values, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(p.columns)))
    ax.set_xticklabels(p.columns, fontsize=10)
    ax.set_yticks(range(len(p)))
    ax.set_yticklabels(p.index, fontsize=8)
    vmax = p.values.max()
    for i in range(len(p)):
        for j in range(len(p.columns)):
            v = int(p.values[i, j])
            ax.text(j, i, v, ha="center", va="center", fontsize=7,
                    color="white" if v < vmax * 0.6 else "black")
    ax.set_title("Significant DE genes per cell type\n(excluded cluster omitted)",
                 fontsize=11)
    fig.colorbar(im, ax=ax, label="n DE genes", fraction=0.04)
    fig.tight_layout()
    out = config.FIG_DIR / "fig4c_de_burden.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.name}")
    return p


def ratio_comparison():
    """Figure: male:female ratio under pooled vs per-cluster analysis."""
    g = pd.read_csv(config.TABLE_DIR / "de_global_summary.csv")
    pooled = dict(zip(g.sex, g.n_sig_guide_threshold))
    pooled_ratio = pooled["Male"] / pooled["Female"]

    s = pd.read_csv(config.TABLE_DIR / "de_per_cluster_summary.csv").dropna(subset=["n_sig"])
    s = s[~s.cluster.astype(str).str.contains("excluded")]
    t = s.groupby("sex")["n_sig"].sum()
    cluster_ratio = t["Male"] / t["Female"]

    labels = ["Pooled\n(this study)", "Per-cluster\n(this study)",
              "Baker et al.\n(reported)"]
    vals = [pooled_ratio, cluster_ratio, 2.15]

    fig, ax = plt.subplots(figsize=(5, 3.8))
    b = ax.bar(labels, vals, color=[RED, BLUE, PURPLE], width=0.6)
    ax.bar_label(b, fmt="%.2f×", fontsize=10, padding=3)
    ax.axhline(1, ls="--", lw=1, c="black")
    ax.text(2.45, 1.04, "no bias", fontsize=7, ha="right")
    ax.set_ylabel("Male : female DE gene ratio", fontsize=9)
    ax.set_ylim(0, max(vals) * 1.35)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Male bias depends on the unit of analysis", fontsize=11)
    fig.tight_layout()
    out = config.FIG_DIR / "fig4d_ratio_comparison.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.name}  (pooled {pooled_ratio:.2f}x, "
          f"per-cluster {cluster_ratio:.2f}x)")
    return pooled, t


def permutation_bars():
    """Figure: treatment effect against its two null contrasts, per sex."""
    p = config.TABLE_DIR / "permutation_control.csv"
    if not p.exists():
        print("  SKIP permutation figure (run scripts/10 first)")
        return
    d = pd.read_csv(p)
    if "axis" not in d.columns:
        print("  SKIP permutation figure (no axis column; re-run scripts/10)")
        return

    ymax = d.n_DE.max() * 1.25
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2), sharey=True)
    for ax, sex in zip(axes, ["Male", "Female"]):
        s = d[d.sex == sex].set_index("axis")["n_DE"]
        order = ["treatment", "replicate", "diagonal"]
        b = ax.bar(range(3), [s.get(o, 0) for o in order],
                   color=[RED, BLUE, GREY], width=0.62)
        ax.bar_label(b, fontsize=10, padding=2)
        ax.set_xticks(range(3))
        ax.set_xticklabels(["treatment\n(real)", "replicate\n(null)",
                            "diagonal\n(null)"], fontsize=9)
        ax.set_ylim(0, ymax)
        ax.set_title(sex, fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
    axes[0].set_ylabel("DE genes", fontsize=9)
    fig.suptitle("Treatment effect versus its own background\n"
                 "Two independent null contrasts per sex, same cells and thresholds",
                 fontsize=11)
    fig.tight_layout()
    out = config.FIG_DIR / "fig7_permutation_control.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.name}")


def volcanoes():
    """Figures: per-sex volcano plots with trimmed axes and stacked labels."""
    for sex in ["male", "female"]:
        f = config.TABLE_DIR / f"de_global_{sex}_all.csv"
        if not f.exists():
            print(f"  SKIP volcano {sex} (table missing)")
            continue
        d = pd.read_csv(f).dropna(subset=["logfoldchanges", "pvals_adj"])
        d["y"] = -np.log10(d.pvals_adj.clip(lower=1e-300))
        sig = (d.pvals_adj < config.PADJ_THRESHOLD) & \
              (d.logfoldchanges.abs() > config.LOG2FC_THRESHOLD)
        up, dn = sig & (d.logfoldchanges > 0), sig & (d.logfoldchanges < 0)

        fig, ax = plt.subplots(figsize=(7, 5.5))
        ax.scatter(d.logfoldchanges[~sig], d.y[~sig], s=4, c="lightgrey",
                   rasterized=True)
        ax.scatter(d.logfoldchanges[up], d.y[up], s=14, c=RED,
                   label=f"up ({int(up.sum())})")
        ax.scatter(d.logfoldchanges[dn], d.y[dn], s=14, c=BLUE,
                   label=f"down ({int(dn.sum())})")
        for i, (_, r) in enumerate(d[sig].nlargest(10, "y").iterrows()):
            ax.annotate(r["names"], (r.logfoldchanges, r.y), fontsize=7,
                        xytext=(9, -4 - i * 13), textcoords="offset points",
                        arrowprops=dict(arrowstyle="-", lw=.4, color="grey"))
        lim = np.ceil(d.logfoldchanges[sig].abs().max()) + 1.5 if sig.any() else 5
        ax.set_xlim(-lim, lim)
        for x in (config.LOG2FC_THRESHOLD, -config.LOG2FC_THRESHOLD):
            ax.axvline(x, ls="--", lw=.7, c="k")
        ax.axhline(-np.log10(config.PADJ_THRESHOLD), ls="--", lw=.7, c="k")
        ax.set_xlabel("log2 fold change (Cocaine / Sucrose)", fontsize=9)
        ax.set_ylabel("-log10 adjusted p", fontsize=9)
        ax.set_title(f"{sex.capitalize()}: cocaine vs. sucrose", fontsize=11)
        ax.legend(frameon=False, fontsize=9, loc="upper left")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        out = config.FIG_DIR / f"05_volcano_{sex}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {out.name}")


def remove_superseded():
    """Delete the scanpy defaults that scripts 02 and 04 rewrite each run."""
    stale = [
        "violin_02_qc_before.png", "violin_02_qc_after.png",
        "scatter_02_counts_vs_genes.png",
        "umap_04_umap_celltypes.png",
        "dotplot__04_canonical_markers_dotplot.png",
        "matrixplot__04_canonical_markers_matrix.png",
        "pca_variance_ratio_03_pca_variance.png",
        "umap_03_umap_clusters.png", "umap_03_umap_metadata.png",
        "05_de_burden_heatmap.png", "05_sex_concordance.png",
    ]
    n = 0
    for f in stale:
        p = config.FIG_DIR / f
        if p.exists():
            p.unlink()
            n += 1
    print(f"  removed {n} superseded default plots")


def main():
    print("1. Removing superseded scanpy defaults")
    remove_superseded()

    print("\n2. QC figure")
    run("make_qc_figure.py")

    print("\n3. Clustering figures")
    run("make_cluster_figure.py")

    print("\n4. Supplementary figures")
    run("make_supp_figures.py")

    print("\n5. Annotation figures")
    run("make_annotation_figure.py")

    print("\n6. Differential expression figures")
    volcanoes()
    p = de_burden_heatmap()
    pooled, per_cluster = ratio_comparison()
    permutation_bars()

    print("\n7. Graphical abstract and flowchart")
    print("   NOTE: numbers in these two are hard-coded. Edit")
    print("   tools/make_graphical_abstract.py with the values below, then run it.")

    print("\n" + "=" * 62)
    print("NUMBERS THAT CHANGED — update the report and the graphical abstract")
    print("=" * 62)
    print(f"  pooled DE:       male {pooled['Male']}, female {pooled['Female']}")
    print(f"  pooled ratio:    {pooled['Male']/pooled['Female']:.2f}x")
    print(f"  per-cluster:     male {per_cluster['Male']}, female {per_cluster['Female']}"
          f"  ({per_cluster['Male']/per_cluster['Female']:.2f}x)")
    perm = config.TABLE_DIR / "permutation_control.csv"
    if perm.exists():
        d = pd.read_csv(perm)
        for sex in ["Male", "Female"]:
            s = d[d.sex == sex].set_index("axis")["n_DE"] if "axis" in d else None
            if s is not None:
                print(f"  permutation {sex:<7} treatment {s.get('treatment')}, "
                      f"replicate {s.get('replicate')}, diagonal {s.get('diagonal')}")
    print("\n  top responding cell types (summed across sexes):")
    for ct, v in p.sum(axis=1).sort_values(ascending=False).head(5).items():
        print(f"    {ct[:44]:<44} {int(v)}")


if __name__ == "__main__":
    main()
