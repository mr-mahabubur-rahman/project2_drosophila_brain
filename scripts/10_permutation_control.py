#!/usr/bin/env python3
"""
10_permutation_control.py — is the differential expression signal above its own
background?

THE PROBLEM THIS ADDRESSES
Reporting "90 differentially expressed genes in males, 14 in females" is only
meaningful against a baseline. With two biological replicates per arm, some
genes will differ between any two samples for reasons unrelated to treatment:
dissection session, dissociation efficiency, sequencing depth, ambient RNA.
Without a null, there is no way to know whether 90 is large or 14 is small.

THE TEST
Within each sex there are four samples: two cocaine, two sucrose. There are
six ways to split them 2-versus-2. One is the true treatment split, one is its
exact complement, and four cut ACROSS the treatment assignment — each pairing
one cocaine sample with one sucrose sample on either side.

Those four are nulls. No real treatment difference exists between the groups
they define, so any genes they return measure between-sample variation at the
same cell numbers, the same test and the same thresholds as the real comparison.

INTERPRETATION
  real >> nulls   the signal is specific to treatment
  real ~ nulls    the signal is indistinguishable from between-sample variation
  real << nulls   between-sample variation exceeds the treatment effect, and the
                  reported gene count cannot be attributed to treatment

This is the strongest permutation test the design permits. With only four
samples per sex there are four null groupings, so the resolution is coarse --
it cannot produce a p-value, only an order-of-magnitude comparison. That is
still far more informative than reporting the gene count alone.

Run:  python scripts/10_permutation_control.py
"""

import itertools
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
    DE_METHOD, LOG2FC_THRESHOLD, PADJ_THRESHOLD,
)

sc.settings.verbosity = 0


def de_count(adata, assignment):
    """DE genes under an arbitrary two-group split of the samples."""
    sub = adata.copy()
    sub.obs["_grp"] = pd.Categorical(sub.obs["sample"].map(assignment))
    if sub.obs["_grp"].nunique() < 2:
        return None, None
    sc.tl.rank_genes_groups(sub, groupby="_grp", reference="B",
                            method=DE_METHOD, use_raw=True)
    d = sc.get.rank_genes_groups_df(sub, group="A").dropna(
        subset=["logfoldchanges", "pvals_adj"])
    sig = d[(d.logfoldchanges.abs() > LOG2FC_THRESHOLD) &
            (d.pvals_adj < PADJ_THRESHOLD)]
    del sub
    return int(len(sig)), set(sig["names"])


def main():
    adata = sc.read_h5ad(H5AD_ANNOTATED)
    print(f"Loaded {adata.n_obs:,} cells\n")

    rows = []
    gene_sets = {}

    for sex in ["Male", "Female"]:
        sub = adata[adata.obs.sex == sex].copy()
        samples = sorted(sub.obs["sample"].unique())
        truth = {s: sub.obs.loc[sub.obs["sample"] == s, "treatment"].iloc[0]
                 for s in samples}

        print("=" * 66)
        print(f"{sex}  ({sub.n_obs:,} cells, {len(samples)} samples)")
        print("=" * 66)

        for combo in itertools.combinations(samples, 2):
            assign = {s: ("A" if s in combo else "B") for s in samples}
            # Is this the true treatment split (or its complement)?
            as_truth = all((assign[s] == "A") == (truth[s] == "Cocaine")
                           for s in samples)
            as_compl = all((assign[s] == "A") == (truth[s] != "Cocaine")
                           for s in samples)
            if as_compl:
                continue                      # same contrast, reversed sign
            # Null splits also come in complementary pairs: {A,B} vs {C,D} is
            # the same contrast as {C,D} vs {A,B}. Keep only one of each, or
            # the figure implies four independent nulls where there are two.
            if not as_truth and samples[0] not in combo:
                continue
            kind = "REAL (cocaine vs sucrose)" if as_truth else "null"

            n, genes = de_count(sub, assign)
            if n is None:
                continue
            label = " + ".join(c.replace(f"{sex}_", "") for c in combo)
            # Which axis does this split follow? The replicate axis groups both
            # R1 samples against both R2 samples.
            reps = {c: c.rsplit("_", 1)[1] for c in samples}
            axis = ("treatment" if as_truth
                    else "replicate" if len({reps[c] for c in combo}) == 1
                    else "diagonal")
            rows.append({"sex": sex, "kind": ("real" if as_truth else "null"), "axis": axis,
                         "group_A": label, "n_DE": n})
            gene_sets[(sex, kind, label)] = genes
            print(f"  {kind:<26} A = {label:<30} {n:>4} genes")

        del sub

    df = pd.DataFrame(rows)
    df.to_csv(TABLE_DIR / "permutation_control.csv", index=False)

    # ---- Summary ---------------------------------------------------------
    print("\n" + "=" * 66)
    print("SUMMARY")
    print("=" * 66)
    for sex in ["Male", "Female"]:
        d = df[df.sex == sex]
        real = d[d.kind == "real"].n_DE.iloc[0]
        nulls = d[d.kind == "null"].n_DE.values
        print(f"\n{sex}:")
        print(f"  real treatment split : {real}")
        print(f"  null splits          : {sorted(nulls)} "
              f"(median {np.median(nulls):.0f})")
        if np.median(nulls) > 0:
            print(f"  signal-to-background : {real/np.median(nulls):.2f}x")
        if real > max(nulls):
            print("  -> signal exceeds every null: specific to treatment")
        elif real < min(nulls):
            print("  -> signal BELOW every null: between-sample variation exceeds")
            print("     the treatment effect; the gene count cannot be attributed")
            print("     to treatment at this threshold")
        else:
            print("  -> signal within the null range: not distinguishable from")
            print("     between-sample variation")

    # ---- Figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=False)
    for ax, sex in zip(axes, ["Male", "Female"]):
        d = df[df.sex == sex]
        real = d[d.kind == "real"].n_DE.iloc[0]
        nulls = d[d.kind == "null"].n_DE.values
        ax.scatter(np.full(len(nulls), 0.6) + np.linspace(-.08, .08, len(nulls)),
                   nulls, s=45, c="grey", zorder=3, label="null splits")
        ax.scatter([1.4], [real], s=110, c="#C44E52", marker="D", zorder=3,
                   label="true treatment split")
        ax.hlines(np.median(nulls), 0.35, 0.85, color="grey", lw=1.5, ls="--")
        ax.set_xlim(0.2, 1.8)
        ax.set_xticks([0.6, 1.4])
        ax.set_xticklabels(["null", "real"], fontsize=9)
        ax.set_ylim(0, max(max(nulls), real) * 1.25)
        ax.set_ylabel("DE genes", fontsize=9)
        ax.set_title(f"{sex}: real = {real}, null median = {np.median(nulls):.0f}",
                     fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("Differential expression against its own background\n"
                 "(four alternative 2-vs-2 sample splits per sex)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig7_permutation_control.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote fig7_permutation_control.png")

    # ---- How much of the real signal is shared with the nulls? -----------
    print("\n" + "=" * 66)
    print("Overlap between the real gene list and the null gene lists")
    print("=" * 66)
    print("Genes appearing in both suggest between-sample variation rather than")
    print("a treatment effect.\n")
    for sex in ["Male", "Female"]:
        real_key = [k for k in gene_sets
                    if k[0] == sex and k[1].startswith("REAL")]
        if not real_key:
            continue
        real_genes = gene_sets[real_key[0]]
        null_union = set().union(*[gene_sets[k] for k in gene_sets
                                   if k[0] == sex and k[1] == "null"])
        shared = real_genes & null_union
        if real_genes:
            print(f"  {sex}: {len(shared)}/{len(real_genes)} real genes also appear "
                  f"in a null split ({100*len(shared)/len(real_genes):.0f}%)")
            if shared:
                print(f"     e.g. {sorted(shared)[:8]}")

    print("\nReport the real count alongside the null range. A gene count given")
    print("without its background is uninterpretable at this replicate number.")


if __name__ == "__main__":
    main()
