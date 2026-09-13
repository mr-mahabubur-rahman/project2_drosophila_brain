"""
01_load_data.py — read the eight 10x matrices and merge them into one AnnData.

WHAT THIS STEP IS FOR
Each sample is a separate CellRanger output directory. Merging them early
(rather than analysing each separately) is what lets us cluster all cells in
one shared embedding, which is the only way to ask "does cluster 17 respond
differently in males and females?".

TWO DATA PROBLEMS HANDLED HERE
Both were found by debugging, not by reading documentation. Both are silent
failures if left alone.

1. SCRAMBLED SAMPLE LABELS. The supplied folder names do not match their
   contents. Sex markers separate the samples by the TREATMENT field of the
   folder name, not the SEX field. Corrected via config.SAMPLE_REMAP, built
   from the deposition order in the authors' published R code. Original folder
   names are preserved in obs['sample_folder'] so the correction is auditable.
   Full reasoning: docs/data_provenance.md and the comment block in config.py.

2. A GENE NAMED 'nan'. Drosophila gene nanchung (Dmel_CG5842) has the symbol
   'nan', which pandas parses as a missing value. The gene ends up with a
   NaN name, which breaks h5ad writing and would silently fail every marker
   lookup downstream. Recovered via config.GENE_NAME_NA_FIXES.

A third problem, space-delimited features.tsv.gz, is fixed upstream by
tools/fix_features_separator.py -- run that once before this script.

Run:  python scripts/01_load_data.py
"""

import sys

import anndata as ad
import pandas as pd
import scanpy as sc

from config import (
    DATA_DIR, H5AD_RAW, FIG_DIR, RANDOM_SEED,
    SAMPLE_REMAP, APPLY_SAMPLE_REMAP, GENE_NAME_NA_FIXES,
)

sc.settings.verbosity = 3
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)


def find_sample_dirs():
    """Locate the 8 sample folders, each holding the 10x triplet."""
    if not DATA_DIR.exists():
        sys.exit(f"ERROR: {DATA_DIR} does not exist.")

    dirs = sorted(
        p for p in DATA_DIR.iterdir()
        if p.is_dir() and p.name != "original" and (
            (p / "matrix.mtx.gz").exists() or (p / "matrix.mtx").exists()
        )
    )
    if not dirs:
        sys.exit(
            f"ERROR: no 10x sample folders found under {DATA_DIR}.\n"
            "Each folder needs barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz."
        )
    return dirs


def sample_metadata(folder_name):
    """
    Resolve a folder name to its TRUE sex, treatment and replicate.

    We deliberately do NOT parse the folder name for sex and treatment -- that
    is exactly what produced the wrong answer. The remap table is the authority.
    """
    if APPLY_SAMPLE_REMAP:
        if folder_name not in SAMPLE_REMAP:
            sys.exit(
                f"ERROR: folder '{folder_name}' is not in config.SAMPLE_REMAP.\n"
                "The remap was built for a specific set of eight folder names.\n"
                "Do not guess -- verify the sample identity before adding it."
            )
        return SAMPLE_REMAP[folder_name]

    # Fallback: trust the folder name. Only for reproducing the original
    # (incorrect) labelling as a comparison.
    parts = folder_name.split("_")
    if len(parts) != 3:
        sys.exit(f"ERROR: '{folder_name}' does not match Sex_Treatment_Replicate.")
    return parts[0], parts[1], parts[2], "unknown"


def fix_na_gene_names(adata):
    """
    Recover gene symbols that pandas parsed as missing values.

    Matching is by gene_id where available, falling back to positional repair,
    because after concatenation the gene_ids column may not survive.
    """
    bad = pd.isna(adata.var_names)
    n_bad = int(bad.sum())
    if n_bad == 0:
        return adata

    print(f"  {n_bad} gene symbol(s) parsed as missing by pandas -- recovering")

    names = pd.Series(adata.var_names, dtype=object)

    if "gene_ids" in adata.var.columns:
        ids = adata.var.loc[bad, "gene_ids"]
        replacements = [
            GENE_NAME_NA_FIXES.get(gid, str(gid).replace("Dmel_", ""))
            for gid in ids
        ]
    else:
        # No gene_ids to match on: use the configured replacements in order.
        replacements = list(GENE_NAME_NA_FIXES.values())[:n_bad]
        if len(replacements) < n_bad:
            sys.exit(
                f"ERROR: {n_bad} unnamed genes but only {len(replacements)} "
                "entries in config.GENE_NAME_NA_FIXES. Investigate before "
                "proceeding -- do not let a gene through unnamed."
            )

    names[bad] = replacements
    adata.var_names = pd.Index(names.astype(str))
    adata.var_names_make_unique()
    print(f"  recovered as: {replacements}")
    return adata


def main():
    sc.settings.seed = RANDOM_SEED
    sample_dirs = find_sample_dirs()
    print(f"Found {len(sample_dirs)} samples")
    if len(sample_dirs) != 8:
        print(f"WARNING: expected 8 samples, found {len(sample_dirs)}.")

    if APPLY_SAMPLE_REMAP:
        print("\nApplying corrected sample labels (config.SAMPLE_REMAP).")
        print("Folder names do NOT match contents -- see docs/data_provenance.md\n")
    else:
        print("\nWARNING: APPLY_SAMPLE_REMAP is False. Using folder names as-is.")
        print("These are known to be incorrect. For comparison only.\n")

    adatas = {}
    for path in sample_dirs:
        folder = path.name
        sex, treatment, replicate, geo = sample_metadata(folder)
        sample_id = f"{sex}_{treatment}_R{replicate}"

        # var_names='gene_symbols' gives readable names (repo, elav, ple)
        # instead of FlyBase IDs. Note this is also what exposes the 'nan'
        # problem -- gene IDs would have avoided it, at the cost of making
        # every marker lookup downstream fail to match.
        a = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
        a.var_names_make_unique()

        a.obs["sample"] = sample_id
        a.obs["sample_folder"] = folder
        a.obs["sex"] = sex
        a.obs["treatment"] = treatment
        a.obs["replicate"] = replicate
        a.obs["geo_position"] = geo
        a.obs["condition"] = f"{sex}_{treatment}"

        arrow = f"{folder:>18} -> {sample_id:<22}" if folder != sample_id else folder
        print(f"  {arrow} ({geo})  {a.n_obs:,} cells x {a.n_vars:,} genes")
        adatas[sample_id] = a

    # join='outer' keeps every gene seen in any sample, filling absent genes
    # with zeros. join='inner' would silently drop genes missing from a single
    # sample -- including, potentially, a marker needed later.
    # index_unique='-' appends the sample key to barcodes so identical
    # barcodes from different runs do not collide.
    adata = ad.concat(adatas, label="sample_batch", index_unique="-", join="outer")
    adata.var_names_make_unique()

    adata = fix_na_gene_names(adata)

    for col in ["sample", "sample_folder", "sex", "treatment", "replicate",
                "geo_position", "condition"]:
        adata.obs[col] = adata.obs[col].astype("category")

    adata.uns["label_correction"] = (
        "Folder names mismatched contents. Corrected from the deposition order "
        "in the authors' published R code; sex confirmed by roX1/roX2 and "
        "Yp1-Yp3/Sxl. Original folder names in obs['sample_folder']."
        if APPLY_SAMPLE_REMAP else "NOT APPLIED -- folder labels used as supplied."
    )

    print(f"\nMerged dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print("(Paper reports 86,224 cells AFTER their QC -- expect more here.)")

    # ---- Sanity checks --------------------------------------------------
    assert not pd.isna(adata.var_names).any(), "unnamed gene survived -- stop"

    print("\nCells per sample:")
    print(adata.obs["sample"].value_counts().sort_index().to_string())

    print("\nDesign (must be a full 2x2x2 -- eight non-zero cells):")
    print(pd.crosstab(adata.obs["sex"],
                      [adata.obs["treatment"], adata.obs["replicate"]]).to_string())

    adata.write(H5AD_RAW)
    print(f"\nWrote {H5AD_RAW}")
    print("\nNext: python tools/verify_remap.py  (confirms the labels hold up)")


if __name__ == "__main__":
    main()
