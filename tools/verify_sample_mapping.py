#!/usr/bin/env python3
"""
verify_sample_mapping.py — confirm your GSM -> sample mapping is not inverted.

WHY THIS MATTERS MORE THAN IT LOOKS
A swapped sex or treatment label produces a pipeline that runs flawlessly,
generates beautiful figures, and answers the wrong question. Nothing crashes.
Nothing looks odd. The error surfaces only when a reviewer notices your
"female-specific" gene is a known male-specific lncRNA.

THE TEST
roX1 and roX2 are lncRNAs of the dosage-compensation complex, expressed almost
exclusively in males. Female samples should show near-zero expression.
If your "Female" samples are roX-high, the sex labels are swapped.

Treatment labels cannot be checked this way -- there is no marker that says
"this fly ate cocaine". For those you must verify against the GEO series page
and Supplemental Table S2 of the paper. This script reminds you to.

Run AFTER step 01:  uv run tools/verify_sample_mapping.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import numpy as np
import pandas as pd
import scanpy as sc

from config import H5AD_RAW, TABLE_DIR, SEX_CONTROL_GENES

sc.settings.verbosity = 0


def main():
    if not H5AD_RAW.exists():
        sys.exit(f"{H5AD_RAW} not found. Run scripts/01_load_data.py first.")

    adata = sc.read_h5ad(H5AD_RAW)
    print(f"Loaded {adata.n_obs:,} cells from {adata.obs['sample'].nunique()} samples\n")

    present = [g for g in SEX_CONTROL_GENES if g in adata.var_names]
    if not present:
        print("WARNING: none of the roX genes were found under any expected name.")
        print(f"Tried: {SEX_CONTROL_GENES}")
        candidates = [v for v in adata.var_names if "roX" in v or "rox" in v.lower()]
        print(f"Genes containing 'roX': {candidates[:10]}")
        print("\nFall back to verifying against the GEO page and Table S2 manually.")
        return

    print(f"Using sex-control genes: {present}\n")

    # Light normalization just for this check -- we do not save it.
    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    df = sc.get.obs_df(a, keys=present)
    df["sample"] = a.obs["sample"].values
    df["sex_label"] = a.obs["sex"].values

    by_sample = df.groupby("sample", observed=True)[present].mean().round(3)
    by_sample["declared_sex"] = (
        df.groupby("sample", observed=True)["sex_label"].first()
    )
    by_sample["roX_mean"] = by_sample[present].mean(axis=1).round(3)

    print("Mean log-normalized expression per sample:")
    print(by_sample.sort_values("declared_sex").to_string())

    by_sample.to_csv(TABLE_DIR / "sample_mapping_sex_check.csv")

    male = by_sample.loc[by_sample["declared_sex"] == "Male", "roX_mean"]
    female = by_sample.loc[by_sample["declared_sex"] == "Female", "roX_mean"]

    print("\n" + "=" * 66)
    if len(male) == 0 or len(female) == 0:
        print("Could not compare — one sex has no samples. Check tools/sample_map.tsv.")
    elif male.mean() > female.mean() * 2:
        print("PASS — samples labelled Male show substantially higher roX")
        print(f"       (male mean {male.mean():.2f} vs female {female.mean():.2f}).")
        print("       Sex labels are consistent with the biology.")
    elif female.mean() > male.mean() * 2:
        print("FAIL — your 'Female' samples are roX-HIGH.")
        print(f"       (female mean {female.mean():.2f} vs male {male.mean():.2f}).")
        print("       The sex labels are almost certainly INVERTED.")
        print("       Fix tools/sample_map.tsv, re-run organize_data.sh and step 01.")
    else:
        print("UNCLEAR — roX levels are similar between the two labelled groups")
        print(f"       (male {male.mean():.2f} vs female {female.mean():.2f}).")
        print("       This is itself suspicious. Inspect the per-sample table above:")
        print("       individual samples should separate cleanly into high and low.")
    print("=" * 66)

    print("\nSTILL TO VERIFY BY HAND (no marker gene can do this for you):")
    print("  Cocaine vs Sucrose assignment, against")
    print("    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152495")
    print("    and Supplemental Table S2 of docs/paper.pdf")
    print("  Record the verification in reports/AI_USAGE_DISCLOSURE.md.")


if __name__ == "__main__":
    main()
