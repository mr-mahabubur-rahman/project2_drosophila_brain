#!/usr/bin/env python3
"""
revert_remap.py — restore the sample labels as deposited.

WHY THIS REVERSION
The sample labels were corrected earlier in this project on the basis of
sex-specific marker genes, using a reconstruction from the deposition order in
the authors' published R code. Supplemental Table S2 of Baker et al. (2021)
subsequently showed that reconstruction to be wrong.

Table S2 lists per-sample cell counts. All eight match the folders exactly:

    Table S2 label      cells     folder as supplied
    -------------------------------------------------------
    F Sucrose Rep 1     9,072     Female_Sucrose_1
    F Sucrose Rep 2    11,693     Female_Sucrose_2
    M Sucrose Rep 1    13,193     Male_Sucrose_1
    M Sucrose Rep 2    11,033     Male_Sucrose_2
    F Cocaine Rep 1    13,072     Female_Cocaine_1
    F Cocaine Rep 2     9,367     Female_Cocaine_2
    M Cocaine Rep 1    10,437     Male_Cocaine_1
    M Cocaine Rep 2    11,124     Male_Cocaine_2
    -------------------------------------------------------
    Total              88,991     (matches merged dataset)

Cell count is intrinsic to each matrix file, so this ties every folder to the
sample name the authors themselves used. It is stronger evidence about labelling
than the deposition-order inference it replaces.

WHAT THIS DOES NOT RESOLVE
The marker contradiction is unchanged and is NOT explained by reverting:

    metric              4 Sucrose samples    4 Cocaine samples
    -----------------------------------------------------------
    roX1 (male)         3.8 - 4.0            0.06 - 0.11
    roX2 (male)         2.1 - 2.4            0.008 - 0.015
    Yp1-Yp3 (female)    0.000 - 0.010        0.028 - 0.422
    Sxl (female)        0.70 - 0.98          1.23 - 1.76

Three independent sex markers separate the samples by TREATMENT, crossing both
declared sexes. Yolk proteins are not transcribed in males, so a value of 0.251
in a sample labelled male is not explicable as regulation.

This is recorded as an unresolved observation, not corrected. Resolving it
requires information the deposited data does not contain and should be raised
with the original authors.

Run:  python tools/revert_remap.py          # show what changes
      python tools/revert_remap.py --apply  # write it
"""

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "scripts" / "config.py"

IDENTITY_REMAP = '''# SAMPLE LABELS — as deposited, verified against Supplemental Table S2
#
# An earlier version of this file remapped these labels on the basis of
# sex-specific markers, using a reconstruction from the deposition order in the
# authors' published R code. Supplemental Table S2 shows that reconstruction was
# wrong: it lists per-sample cell counts, and all eight match the folders as
# supplied (9,072 / 11,693 / 13,193 / 11,033 / 13,072 / 9,367 / 10,437 / 11,124;
# total 88,991). Cell count is intrinsic to each matrix, so this ties every
# folder to the name the authors used. The labels below are therefore the
# deposited ones, unchanged.
#
# UNRESOLVED: sex-specific markers do not agree with these labels. roX1/roX2
# (male-specific) read 3.8-4.0 in all four Sucrose samples and 0.06-0.11 in all
# four Cocaine samples; Yp1-Yp3 (female-specific) do the reverse. The split
# follows treatment, crossing both declared sexes. Yolk proteins are not
# transcribed in males, so this is not explicable as regulation. Documented in
# docs/data_provenance.md; not corrected here, because correcting it would
# require information the deposited data does not contain.
#
# folder name -> (sex, treatment, replicate, GEO position)
SAMPLE_REMAP = {
    "Female_Sucrose_1": ("Female", "Sucrose", "1", "S1"),
    "Female_Sucrose_2": ("Female", "Sucrose", "2", "S2"),
    "Male_Sucrose_1":   ("Male",   "Sucrose", "1", "S3"),
    "Male_Sucrose_2":   ("Male",   "Sucrose", "2", "S4"),
    "Female_Cocaine_1": ("Female", "Cocaine", "1", "S5"),
    "Female_Cocaine_2": ("Female", "Cocaine", "2", "S6"),
    "Male_Cocaine_1":   ("Male",   "Cocaine", "1", "S7"),
    "Male_Cocaine_2":   ("Male",   "Cocaine", "2", "S8"),
}

# True = use the table above (labels as deposited). The table is now an identity
# mapping, so this switch no longer changes any assignment; it is retained so
# 01_load_data.py continues to work unmodified.
APPLY_SAMPLE_REMAP = True
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    src = CONFIG.read_text(encoding="utf-8")

    start = src.find("# SAMPLE LABEL CORRECTION")
    if start == -1:
        start = src.find("# folder name -> (sex, treatment, replicate, GEO position)")
    end = src.find("APPLY_SAMPLE_REMAP")
    if start == -1 or end == -1:
        sys.exit("Could not locate the SAMPLE_REMAP block in config.py. "
                 "Edit it by hand using the table in this script's docstring.")
    end = src.find("\n", src.find("\n", end) + 1)

    print("REPLACING this block in scripts/config.py:\n")
    print("-" * 66)
    old = src[start:end]
    print(old[:700] + ("\n  ...\n" if len(old) > 700 else ""))
    print("-" * 66)
    print("\nWITH an identity mapping (labels as deposited).\n")

    print("Mapping after reversion:")
    for folder, sex, trt in [
        ("Female_Sucrose_1", "Female", "Sucrose"), ("Female_Sucrose_2", "Female", "Sucrose"),
        ("Male_Sucrose_1", "Male", "Sucrose"), ("Male_Sucrose_2", "Male", "Sucrose"),
        ("Female_Cocaine_1", "Female", "Cocaine"), ("Female_Cocaine_2", "Female", "Cocaine"),
        ("Male_Cocaine_1", "Male", "Cocaine"), ("Male_Cocaine_2", "Male", "Cocaine"),
    ]:
        print(f"  {folder:<20} -> {sex:<7} {trt}")

    if not args.apply:
        print("\nReport only. Re-run with --apply to write the change.")
        return

    backup = CONFIG.with_suffix(".py.pre-revert")
    shutil.copy2(CONFIG, backup)
    CONFIG.write_text(src[:start] + IDENTITY_REMAP + src[end:], encoding="utf-8")
    print(f"\nBacked up to {backup.name}")
    print(f"Rewrote {CONFIG}")

    print("""
NEXT — the whole pipeline must be re-run. Every checkpoint downstream of
step 01 carries the old labels:

    python scripts/01_load_data.py
    python scripts/02_qc_filter.py
    python scripts/03_normalize_cluster.py
    python scripts/04_annotate_clusters.py
    python scripts/05_de_analysis.py
    python scripts/05b_pseudobulk_check.py
    python scripts/07_interaction_test.py
    python scripts/08_sensitivity.py
    python scripts/09_reproduce_fig3.py
    python scripts/10_permutation_control.py

Clustering and annotation will be unchanged -- neither depends on the labels.
Every differential expression result will change, because the comparison groups
are now different cells.
""")


if __name__ == "__main__":
    main()
