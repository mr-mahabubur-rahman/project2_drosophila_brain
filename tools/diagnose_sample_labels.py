#!/usr/bin/env python3
"""
diagnose_sample_labels.py — determine which metadata field actually encodes sex.

WHY THIS EXISTS
tools/verify_sample_mapping.py found that roX1/roX2 separate the eight samples
cleanly into two groups of four -- but the split follows the TREATMENT label,
not the SEX label. Since roX1/roX2 are dosage-compensation lncRNAs and among the
most male-specific transcripts in Drosophila, that pattern is very hard to
explain as a drug effect and easy to explain as mislabelling.

One marker pair is not enough to act on. This script tests several independent
sex markers, and reports how each groups the samples, so the conclusion rests on
converging evidence rather than a single gene.

MARKERS USED
  Male-biased:
    lncRNA:roX1, lncRNA:roX2   dosage compensation, male-specific
    mle, msl-1, msl-2, msl-3   MSL complex; msl-2 is translated only in males
  Female-biased:
    Sxl                        Sex-lethal; functional protein only in females
    Yp1, Yp2, Yp3              yolk proteins; female fat body, carried over in
                               whole-brain dissections
    tra                        transformer; female-specific splice form

If the male-biased and female-biased panels give OPPOSITE groupings of the same
four samples, the sex assignment is unambiguous.

Run:  python tools/diagnose_sample_labels.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import numpy as np
import pandas as pd
import scanpy as sc

from config import H5AD_RAW, TABLE_DIR

sc.settings.verbosity = 0

MALE_MARKERS = ["lncRNA:roX1", "lncRNA:roX2", "roX1", "roX2",
                "mle", "msl-1", "msl-2", "msl-3"]
FEMALE_MARKERS = ["Sxl", "Yp1", "Yp2", "Yp3", "tra"]


def present(adata, genes):
    return [g for g in genes if g in adata.var_names]


def main():
    if not H5AD_RAW.exists():
        sys.exit(f"{H5AD_RAW} not found. Run notebook/script 01 first.")

    adata = sc.read_h5ad(H5AD_RAW)
    print(f"Loaded {adata.n_obs:,} cells from {adata.obs['sample'].nunique()} samples\n")

    male_g = present(adata, MALE_MARKERS)
    female_g = present(adata, FEMALE_MARKERS)
    print(f"Male-biased markers found:   {male_g}")
    print(f"Female-biased markers found: {female_g}")
    if not male_g or not female_g:
        print("\nWARNING: one panel is empty; the comparison will be weak.")

    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    all_g = male_g + female_g
    df = sc.get.obs_df(a, keys=all_g)
    df["sample"] = a.obs["sample"].values

    tbl = df.groupby("sample", observed=True)[all_g].mean().round(3)

    # Parse the folder name back into its two label fields.
    tbl["field1_sex_label"] = [s.split("_")[0] for s in tbl.index]
    tbl["field2_treatment_label"] = [s.split("_")[1] for s in tbl.index]

    if male_g:
        tbl["MALE_panel"] = tbl[male_g].mean(axis=1).round(3)
    if female_g:
        tbl["FEMALE_panel"] = tbl[female_g].mean(axis=1).round(3)

    print("\n" + "=" * 78)
    print("PER-SAMPLE MARKER EXPRESSION (log-normalized mean)")
    print("=" * 78)
    cols = all_g + [c for c in ["MALE_panel", "FEMALE_panel"] if c in tbl]
    print(tbl[cols].to_string())

    tbl.to_csv(TABLE_DIR / "sample_label_diagnosis.csv")

    # ---- Which label field do the markers agree with? -------------------
    print("\n" + "=" * 78)
    print("WHICH LABEL FIELD DOES EACH PANEL FOLLOW?")
    print("=" * 78)

    for panel in ["MALE_panel", "FEMALE_panel"]:
        if panel not in tbl:
            continue
        print(f"\n{panel}:")
        for field in ["field1_sex_label", "field2_treatment_label"]:
            grouped = tbl.groupby(field)[panel].agg(["mean", "min", "max"]).round(3)
            print(f"  grouped by {field}:")
            print(grouped.to_string().replace("\n", "\n    "))

            # Separation ratio: how cleanly do the two groups split?
            means = grouped["mean"]
            if len(means) == 2 and means.min() > 0:
                ratio = means.max() / means.min()
                # A clean split = the groups do not overlap at all
                lo_label, hi_label = means.idxmin(), means.idxmax()
                clean = grouped.loc[hi_label, "min"] > grouped.loc[lo_label, "max"]
                print(f"    ratio {ratio:.1f}x, "
                      f"{'NON-OVERLAPPING' if clean else 'overlapping'}")

    # ---- Verdict ---------------------------------------------------------
    print("\n" + "=" * 78)
    print("READING THIS OUTPUT")
    print("=" * 78)
    print("""
The panel that separates samples into two NON-OVERLAPPING groups of four,
with a large ratio, is tracking real biology. The label field it groups by
is the field that actually encodes sex.

  - Both panels follow field1 (the sex label), in opposite directions
        -> labels are correct; the roX result was a false alarm.

  - Both panels follow field2 (the treatment label), in opposite directions
        -> the two metadata fields are swapped. The label written as
           "treatment" is really sex, and vice versa.

  - Panels disagree, or neither separates cleanly
        -> do not guess. Take it to your instructor.

If the male and female panels both point the same way (e.g. both high in the
same four samples), something is wrong with the markers or the data, not the
labels -- tell me before acting on it.
""")

    print(f"Table written to {TABLE_DIR / 'sample_label_diagnosis.csv'}")
    print("\nDo not relabel anything yet. Confirm with your instructor first --")
    print("they assigned these names, and they can check against the GEO record")
    print("faster than you can infer it.")


if __name__ == "__main__":
    main()
