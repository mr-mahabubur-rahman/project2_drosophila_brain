#!/usr/bin/env python3
"""
08_sensitivity.py — two sensitivity analyses that test the main claims directly.

ANALYSIS A: DOES THE MALE BIAS SURVIVE REMOVING MITOCHONDRIAL GENES?
The central claim is that the pooled male result (90 genes vs 14 in females,
6.43x) is inflated by a mitochondrial signal that differs between treatment
arms in males (median 0.373% vs 0.058%) but not in females. That claim predicts
something specific and falsifiable: removing mitochondrial and respiratory-chain
genes should shrink the male count disproportionately.

If the ratio collapses toward the per-cluster value (1.65x) or the published
value (2.15x), the argument is demonstrated rather than merely asserted. If it
barely moves, the mitochondrial genes are a symptom rather than the cause, and
the claim in the Discussion must be weakened accordingly. Either outcome is
worth reporting.

ANALYSIS B: DOES THE RESPONSE RANKING SURVIVE EQUAL CLUSTER SIZE?
Cluster size correlates with DE gene count (Spearman rho = 0.581), so a raw
ranking partly measures statistical power rather than biological response.
Per-1000-cell normalisation is one correction, but it assumes DE count scales
linearly with cell number, which it does not -- power grows sublinearly.

The cleaner approach is to downsample every cluster to the same number of cells
and re-test. Any remaining differences reflect response magnitude at matched
power. Repeated over several random seeds to show the ranking is stable rather
than an artefact of one draw.

Run:  python scripts/08_sensitivity.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from config import (
    H5AD_ANNOTATED, TABLE_DIR, FIG_DIR,
    DE_METHOD, LOG2FC_THRESHOLD, PADJ_THRESHOLD, RANDOM_SEED,
)

sc.settings.verbosity = 0

N_DOWNSAMPLE = 400      # cells per cluster per analysis; raise if you have time
N_SEEDS = 3             # repeats, to show the ranking is not one lucky draw
MIN_PER_ARM = 50        # skip clusters too small to downsample meaningfully


def de_count(sub):
    """Number of significant cocaine-vs-sucrose genes in a subset."""
    counts = sub.obs["treatment"].value_counts()
    if min(counts.get("Cocaine", 0), counts.get("Sucrose", 0)) < 10:
        return None
    try:
        sc.tl.rank_genes_groups(sub, groupby="treatment", reference="Sucrose",
                                method=DE_METHOD, use_raw=True)
        d = sc.get.rank_genes_groups_df(sub, group="Cocaine").dropna(
            subset=["logfoldchanges", "pvals_adj"])
    except Exception:
        return None
    return int(((d.logfoldchanges.abs() > LOG2FC_THRESHOLD) &
                (d.pvals_adj < PADJ_THRESHOLD)).sum())


# ======================================================================
# A. Mitochondrial / ribosomal exclusion
# ======================================================================
def analysis_a(adata):
    print("=" * 70)
    print("A. Does the male bias survive removing mitochondrial genes?")
    print("=" * 70)

    names = adata.raw.var_names
    mito = names.str.startswith("mt:")
    # Nuclear-encoded respiratory chain and mitochondrial-import genes that
    # rose together in the male result.
    resp_prefixes = ("UQCR", "COX", "ND", "ATPsyn", "Cyt-c", "Mic", "sun",
                     "SdhA", "SdhB", "SdhC", "SdhD", "Vha")
    resp = names.str.startswith(resp_prefixes)
    ribo = names.str.startswith(("RpS", "RpL", "mRpS", "mRpL"))

    sets = {
        "all genes": np.zeros(len(names), dtype=bool),
        "minus mito": mito,
        "minus mito+respiratory": mito | resp,
        "minus mito+respiratory+ribosomal": mito | resp | ribo,
    }

    rows = []
    for label, drop in sets.items():
        keep = ~drop
        dropped_names = set(names[drop])
        print(f"\n{label}  ({int(keep.sum()):,} genes, "
              f"{len(dropped_names)} excluded)")
        counts = {}
        for sex in ["Male", "Female"]:
            sub = adata[adata.obs.sex == sex].copy()
            # rank_genes_groups reads .raw, so restrict .raw ONCE. Building the
            # mask from sub.raw.var_names each time avoids any mismatch if the
            # gene set has already been altered.
            r = sub.raw.to_adata()
            r = r[:, ~r.var_names.isin(dropped_names)].copy()
            sub.raw = r
            n = de_count(sub)
            counts[sex] = n
            print(f"  {sex}: {n} DE genes")
            del sub, r
        ratio = counts["Male"] / counts["Female"] if counts["Female"] else np.nan
        rows.append({"gene_set": label, "n_genes": int(keep.sum()),
                     "male": counts["Male"], "female": counts["Female"],
                     "ratio": round(ratio, 2) if ratio == ratio else None})
        print(f"  ratio = {ratio:.2f}x" if ratio == ratio else "  ratio undefined")

    df = pd.DataFrame(rows)
    df.to_csv(TABLE_DIR / "sensitivity_gene_exclusion.csv", index=False)
    print("\n" + df.to_string(index=False))
    print("\nInterpretation: if the ratio falls toward the per-cluster value")
    print("(1.65x) or the published value (2.15x), the mitochondrial signal was")
    print("driving the pooled male bias. If it holds near 6.4x, weaken the claim.")

    fig, ax = plt.subplots(figsize=(6.5, 4))
    x = range(len(df))
    b = ax.bar(x, df.ratio.fillna(0), color="#4C72B0", width=0.6)
    ax.bar_label(b, fmt="%.2f×", fontsize=9, padding=3)
    ax.axhline(2.15, ls="--", c="#8172B2", lw=1.2, label="Baker et al. 2.15×")
    ax.axhline(1.65, ls=":", c="#55A868", lw=1.2, label="per-cluster 1.65×")
    ax.set_xticks(list(x))
    ax.set_xticklabels([s.replace("minus ", "−\n") for s in df.gene_set],
                       fontsize=8)
    ax.set_ylabel("Male : female DE gene ratio", fontsize=9)
    ax.set_title("Effect of excluding mitochondrial genes", fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig5b_gene_exclusion.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("\nwrote fig5b_gene_exclusion.png")


# ======================================================================
# B. Equal-size downsampling
# ======================================================================
def analysis_b(adata):
    print("\n" + "=" * 70)
    print(f"B. Cluster response ranking at matched power ({N_DOWNSAMPLE} cells each)")
    print("=" * 70)

    types = [c for c in adata.obs["cell_type"].cat.categories
             if "excluded" not in str(c)]
    results = {}

    for seed in range(N_SEEDS):
        rng = np.random.default_rng(RANDOM_SEED + seed)
        print(f"\nseed {seed}:")
        for ct in types:
            m = (adata.obs["cell_type"] == ct).values
            obs = adata.obs[m]
            coc = obs.index[obs.treatment == "Cocaine"]
            suc = obs.index[obs.treatment == "Sucrose"]
            if len(coc) < MIN_PER_ARM or len(suc) < MIN_PER_ARM:
                continue
            k = min(N_DOWNSAMPLE // 2, len(coc), len(suc))
            pick = np.concatenate([rng.choice(coc, k, replace=False),
                                   rng.choice(suc, k, replace=False)])
            sub = adata[pick].copy()
            n = de_count(sub)
            del sub
            if n is None:
                continue
            results.setdefault(str(ct), []).append(n)
            print(f"  {str(ct)[:44]:<44} {n:>4}  (n={2*k})")

    if not results:
        print("No clusters large enough to downsample.")
        return

    df = pd.DataFrame({
        "cell_type": list(results),
        "mean_DE": [np.mean(v) for v in results.values()],
        "sd_DE": [np.std(v) for v in results.values()],
        "n_seeds": [len(v) for v in results.values()],
    }).sort_values("mean_DE", ascending=False)
    df.to_csv(TABLE_DIR / "sensitivity_downsampled_ranking.csv", index=False)

    print("\nRanking at equal cluster size:")
    print(df.round(1).to_string(index=False))

    fig, ax = plt.subplots(figsize=(7, max(3.5, 0.32 * len(df))))
    y = range(len(df))
    ax.barh(list(y), df.mean_DE, xerr=df.sd_DE, color="#4C72B0",
            error_kw=dict(lw=0.8, capsize=2))
    ax.set_yticks(list(y))
    ax.set_yticklabels([c[:44] for c in df.cell_type], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f"DE genes at {N_DOWNSAMPLE} cells per cluster "
                  f"(mean ± SD, {N_SEEDS} seeds)", fontsize=9)
    ax.set_title("Cocaine response at matched statistical power", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig5c_downsampled_ranking.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("\nwrote fig5c_downsampled_ranking.png")
    print("\nCompare against the raw ranking (Kenyon cells first) and the")
    print("per-1000-cell ranking (surface glia first). Agreement between the")
    print("downsampled and normalised rankings would be strong evidence that")
    print("the ordering reflects response rather than power.")


def main():
    adata = sc.read_h5ad(H5AD_ANNOTATED)
    print(f"Loaded {adata.n_obs:,} cells\n")
    analysis_a(adata)
    analysis_b(adata)


if __name__ == "__main__":
    main()
