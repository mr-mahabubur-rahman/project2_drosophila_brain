#!/usr/bin/env python3
"""
remap_sample_labels.py — apply the corrected sex/treatment labels.

THE BUG
The supplied folder names do not match their contents. Sex-specific markers
(roX1/roX2 male; Yp1-Yp3, Sxl female) separate the eight samples cleanly into
two groups of four, but the split follows the TREATMENT field of the folder
name, not the SEX field.

THE CAUSE
The authors' published analysis code
(github.com/vshanka23/The-Drosophila-Brain-on-Cocaine-at-Single-Cell-Resolution)
lists the deposition order of the eight samples:

    S1 Female_sucrose_R1    S5 Female_cocaine_R1
    S2 Female_sucrose_R2    S6 Female_cocaine_R2
    S3 Male_sucrose_R1      S7 Male_cocaine_R1
    S4 Male_sucrose_R2      S8 Male_cocaine_R2

Sorting the supplied folder names alphabetically and pairing them against that
order reproduces the observed marker pattern exactly, for all eight samples.
The folders were evidently built by assigning GSM accessions in GEO order to
folder names in alphabetical order -- two different orderings.

CONFIDENCE
  Sex assignment:       CONFIRMED by two independent marker panels.
  Treatment assignment: INFERRED from the authors' deposition order. No marker
                        gene can report whether a fly ate cocaine, so this rests
                        on the GEO ordering rather than on the data itself.
                        verify_remap.py tests it against an independent
                        prediction from the paper.

WHAT THIS SCRIPT DOES
Rewrites obs columns in the checkpoint from step 01. The original folder labels
are preserved in `sample_folder`, `label_field1`, `label_field2` so nothing is
destroyed and the correction stays auditable.

Run:  python tools/remap_sample_labels.py           # report only
      python tools/remap_sample_labels.py --apply   # rewrite the checkpoint
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import pandas as pd
import scanpy as sc

from config import H5AD_RAW, TABLE_DIR

sc.settings.verbosity = 0

# folder name -> (true sex, true treatment, true replicate, GEO position)
REMAP = {
    "Female_Cocaine_1": ("Female", "Sucrose", "1", "S1"),
    "Female_Cocaine_2": ("Female", "Sucrose", "2", "S2"),
    "Female_Sucrose_1": ("Male",   "Sucrose", "1", "S3"),
    "Female_Sucrose_2": ("Male",   "Sucrose", "2", "S4"),
    "Male_Cocaine_1":   ("Female", "Cocaine", "1", "S5"),
    "Male_Cocaine_2":   ("Female", "Cocaine", "2", "S6"),
    "Male_Sucrose_1":   ("Male",   "Cocaine", "1", "S7"),
    "Male_Sucrose_2":   ("Male",   "Cocaine", "2", "S8"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    if not H5AD_RAW.exists():
        sys.exit(f"{H5AD_RAW} not found. Run step 01 first.")

    adata = sc.read_h5ad(H5AD_RAW)
    print(f"Loaded {adata.n_obs:,} cells\n")

    folders = set(adata.obs["sample"].astype(str).unique())
    unknown = folders - set(REMAP)
    if unknown:
        sys.exit(f"Folders not in the remap table: {sorted(unknown)}\n"
                 "Do not proceed -- the mapping was built for a specific set of names.")

    rows = []
    for folder, (sex, treat, rep, geo) in REMAP.items():
        n = int((adata.obs["sample"].astype(str) == folder).sum())
        f1, f2, f3 = folder.split("_")
        rows.append({
            "folder": folder, "n_cells": n, "geo_position": geo,
            "label_said_sex": f1, "label_said_treatment": f2,
            "TRUE_sex": sex, "TRUE_treatment": treat, "TRUE_replicate": rep,
            "sex_changed": f1 != sex, "treatment_changed": f2 != treat,
        })

    tbl = pd.DataFrame(rows).sort_values("geo_position")
    print(tbl.to_string(index=False))
    tbl.to_csv(TABLE_DIR / "sample_label_remap.csv", index=False)

    print(f"\nSex labels changed:       {int(tbl['sex_changed'].sum())}/8")
    print(f"Treatment labels changed: {int(tbl['treatment_changed'].sum())}/8")

    print("\nCorrected design:")
    design = tbl.groupby(["TRUE_sex", "TRUE_treatment"], observed=True).agg(
        n_samples=("folder", "size"), n_cells=("n_cells", "sum"))
    print(design.to_string())

    if not args.apply:
        print("\nReport only. Re-run with --apply to rewrite the checkpoint.")
        return

    # ---- Apply ----------------------------------------------------------
    obs = adata.obs
    folder_series = obs["sample"].astype(str)

    # Preserve what was there, so the correction is auditable and reversible.
    obs["sample_folder"] = folder_series.values
    obs["label_field1"] = [s.split("_")[0] for s in folder_series]
    obs["label_field2"] = [s.split("_")[1] for s in folder_series]

    obs["sex"] = [REMAP[s][0] for s in folder_series]
    obs["treatment"] = [REMAP[s][1] for s in folder_series]
    obs["replicate"] = [REMAP[s][2] for s in folder_series]
    obs["geo_position"] = [REMAP[s][3] for s in folder_series]
    obs["condition"] = [f"{REMAP[s][0]}_{REMAP[s][1]}" for s in folder_series]
    # Canonical sample name, with R to distinguish it from the folder name.
    obs["sample"] = [f"{REMAP[s][0]}_{REMAP[s][1]}_R{REMAP[s][2]}"
                     for s in folder_series]

    for col in ["sample", "sex", "treatment", "replicate", "condition",
                "sample_folder", "label_field1", "label_field2", "geo_position"]:
        obs[col] = obs[col].astype("category")

    adata.obs = obs
    adata.uns["label_correction"] = (
        "Folder names mismatched contents; corrected using the deposition order "
        "in the authors' published R code. Sex confirmed by roX1/roX2 and "
        "Yp1-Yp3/Sxl; treatment inferred from GEO ordering. Original folder "
        "labels retained in obs['sample_folder']."
    )

    adata.write(H5AD_RAW)
    print(f"\nRewrote {H5AD_RAW}")
    print("\nCells per corrected condition:")
    print(adata.obs["condition"].value_counts().sort_index().to_string())
    print("\nNext: python tools/verify_remap.py")


if __name__ == "__main__":
    main()
