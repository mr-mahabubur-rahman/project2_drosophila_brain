#!/usr/bin/env python3
"""
verify_remap.py — test the corrected labels against independent predictions.

WHY THIS IS NEEDED
remap_sample_labels.py fixed sex using marker genes (solid) and treatment using
the authors' deposition order (an inference). An inference deserves a test that
did not go into making it.

THREE CHECKS, IN DESCENDING ORDER OF STRENGTH

1. SEX MARKERS -- confirmatory, not independent.
   roX1/roX2 should now be high in samples labelled Male, Yp1-Yp3 high in
   Female. This must pass; it is what the sex correction was built on. Failing
   here means the remap was applied wrongly.

2. MALE-BIASED RESPONSE -- independent of how we assigned treatment.
   The paper reports 691 DE genes in males versus 322 in females (~2.1x).
   Under the corrected labels, males should show more cocaine-responsive genes
   than females. Nothing about this prediction was used to build the mapping,
   so agreement is real evidence. Note it tests the sex/treatment PAIRING, not
   which arm is the drug.

3. NAMED GENES -- weak, directional.
   The paper names Rpl41, IA-2, CR34335 and mt:lrRNA as up after cocaine, and
   roX2 and ninaE as down. If the treatment direction were reversed, these would
   all point the wrong way. Circular in the sense that we are using the paper's
   result to check labels derived partly from the paper's repository -- treat as
   supporting, not decisive.

WHAT FAILURE MEANS
If check 2 comes out backwards (females showing far more DE genes than males),
the sex/treatment pairing is wrong and you should stop. If check 3 is reversed
but check 2 holds, the two treatment arms may be swapped -- the DE gene SETS
would still be right, only the sign of every fold change flipped.

Run:  python tools/verify_remap.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import numpy as np
import pandas as pd
import scanpy as sc

from config import (
    H5AD_RAW, TABLE_DIR, LOG2FC_THRESHOLD, PADJ_THRESHOLD,
    PAPER_CORE_RESPONSE_GENES,
)

sc.settings.verbosity = 0


def main():
    adata = sc.read_h5ad(H5AD_RAW)

    if "sample_folder" not in adata.obs:
        sys.exit("Checkpoint has not been remapped. "
                 "Run: python tools/remap_sample_labels.py --apply")

    print(f"Loaded {adata.n_obs:,} cells with corrected labels\n")
    print(pd.crosstab(adata.obs["sex"], adata.obs["treatment"]).to_string())

    # Light QC so the test is not dominated by empty droplets.
    sc.pp.filter_cells(adata, min_genes=300)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print(f"\nAfter light QC: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ---- CHECK 1: sex markers -------------------------------------------
    print("\n" + "=" * 72)
    print("CHECK 1 — sex markers (confirmatory)")
    print("=" * 72)

    male_g = [g for g in ["lncRNA:roX1", "lncRNA:roX2", "roX1", "roX2"]
              if g in adata.var_names]
    female_g = [g for g in ["Yp1", "Yp2", "Yp3", "Sxl"] if g in adata.var_names]

    df = sc.get.obs_df(adata, keys=male_g + female_g)
    df["sex"] = adata.obs["sex"].values
    by_sex = df.groupby("sex", observed=True).mean().round(3)
    print(by_sex.to_string())

    m_high = by_sex.loc["Male", male_g].mean() > by_sex.loc["Female", male_g].mean()
    f_high = by_sex.loc["Female", female_g].mean() > by_sex.loc["Male", female_g].mean()
    check1 = m_high and f_high
    print(f"\n  male markers higher in Male:     {m_high}")
    print(f"  female markers higher in Female: {f_high}")
    print(f"  -> {'PASS' if check1 else 'FAIL'}")

    # ---- CHECK 2: male-biased response ----------------------------------
    print("\n" + "=" * 72)
    print("CHECK 2 — male-biased cocaine response (independent)")
    print("=" * 72)
    print("Paper: 691 DE genes in males vs 322 in females (~2.1x).\n")

    counts = {}
    de_tables = {}
    for sex in ["Male", "Female"]:
        sub = adata[adata.obs["sex"] == sex].copy()
        sc.tl.rank_genes_groups(sub, groupby="treatment", reference="Sucrose",
                                method="wilcoxon")
        d = sc.get.rank_genes_groups_df(sub, group="Cocaine").dropna(
            subset=["logfoldchanges", "pvals_adj"])
        de_tables[sex] = d
        sig = d[(d["logfoldchanges"].abs() > LOG2FC_THRESHOLD) &
                (d["pvals_adj"] < PADJ_THRESHOLD)]
        counts[sex] = len(sig)
        print(f"  {sex}: {sub.n_obs:,} cells -> {len(sig)} DE genes")
        del sub

    if counts["Female"] > 0:
        ratio = counts["Male"] / counts["Female"]
        print(f"\n  Male:Female ratio = {ratio:.2f}x (paper ~2.1x)")
    else:
        ratio = np.inf
        print("\n  Female count is zero -- ratio undefined")

    check2 = counts["Male"] > counts["Female"]
    print(f"  -> {'PASS' if check2 else 'FAIL'} "
          f"({'males respond more, as reported' if check2 else 'DIRECTION IS BACKWARDS'})")

    # ---- CHECK 3: named genes -------------------------------------------
    print("\n" + "=" * 72)
    print("CHECK 3 — genes the paper names (supporting)")
    print("=" * 72)

    rows = []
    for sex in ["Male", "Female"]:
        d = de_tables[sex].set_index("names")
        for direction in ["up", "down"]:
            for g in PAPER_CORE_RESPONSE_GENES[direction]:
                for candidate in (g, f"lncRNA:{g}", g.replace("lncRNA:", "")):
                    if candidate in d.index:
                        lfc = float(d.loc[candidate, "logfoldchanges"])
                        rows.append({
                            "gene": candidate, "sex": sex,
                            "paper_says": direction,
                            "our_log2FC": round(lfc, 3),
                            "agrees": (lfc > 0) if direction == "up" else (lfc < 0),
                        })
                        break

    if rows:
        chk = pd.DataFrame(rows)
        print(chk.to_string(index=False))
        agree_rate = chk["agrees"].mean()
        print(f"\n  agreement: {agree_rate:.0%} of {len(chk)} observations")
        print("  -> above 50% supports the treatment direction; "
              "far below suggests the arms are swapped")
        chk.to_csv(TABLE_DIR / "remap_verification_genes.csv", index=False)
    else:
        agree_rate = None
        print("  None of the named genes were found.")

    # ---- Verdict ---------------------------------------------------------
    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)
    if check1 and check2:
        print("The corrected labels hold up. Sex is confirmed by markers, and the")
        print("male-biased response reproduces independently.")
        if agree_rate is not None and agree_rate < 0.5:
            print("\nBut check 3 disagrees: the two treatment arms may be swapped.")
            print("DE gene SETS would still be correct; only the sign of every")
            print("fold change would flip. Say so explicitly in your report.")
        print("\nProceed to notebook 02. Record all of this in")
        print("docs/data_provenance.md -- the correction is part of your Methods.")
    elif not check1:
        print("Check 1 failed. The remap was not applied correctly. Stop and")
        print("re-run tools/remap_sample_labels.py --apply.")
    else:
        print("Check 1 passed but check 2 is backwards. The sex/treatment pairing")
        print("is wrong. Do not proceed -- the reconstruction needs rethinking.")


if __name__ == "__main__":
    main()
