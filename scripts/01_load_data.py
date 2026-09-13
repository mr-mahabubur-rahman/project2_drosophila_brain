"""
01_load_data.py — read the eight 10x matrices and merge them into one AnnData.

WHAT THIS STEP IS FOR
Each sample is a separate CellRanger output directory. Merging them early
(rather than analysing each separately) is what lets us cluster all cells in
one shared embedding, which is the only way to ask "does cluster 17 respond
differently in males and females?".

The critical part is the METADATA. If sex/treatment/replicate are not attached
to every cell here, nothing downstream works.

Run:  uv run scripts/01_load_data.py
"""

import sys
import anndata as ad
import scanpy as sc

from config import (
    DATA_DIR, H5AD_RAW, FIG_DIR, RANDOM_SEED,
)

sc.settings.verbosity = 3
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)


def find_sample_dirs():
    """Locate the 8 sample folders, each holding the 10x triplet."""
    if not DATA_DIR.exists():
        sys.exit(f"ERROR: {DATA_DIR} does not exist. Run scripts/00_download_data.sh first.")

    dirs = sorted(
        p for p in DATA_DIR.iterdir()
        if p.is_dir() and (
            (p / "matrix.mtx.gz").exists() or (p / "matrix.mtx").exists()
        )
    )
    if not dirs:
        sys.exit(
            f"ERROR: no 10x sample folders found under {DATA_DIR}.\n"
            "Each folder must contain barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz.\n"
            "Run scripts/00_download_data.sh, then check the folder names."
        )
    return dirs


def parse_metadata(sample_name):
    """
    'Female_Cocaine_1' -> ('Female', 'Cocaine', '1')

    We validate rather than trust. A silently mis-parsed folder name would
    put male cells in the female DE analysis and you would never notice.
    """
    parts = sample_name.split("_")
    if len(parts) != 3:
        sys.exit(
            f"ERROR: folder '{sample_name}' does not match Sex_Treatment_Replicate.\n"
            "Rename it, e.g. Female_Cocaine_1."
        )
    sex, treatment, replicate = parts
    if sex not in {"Female", "Male"}:
        sys.exit(f"ERROR: unexpected sex '{sex}' in '{sample_name}'.")
    if treatment not in {"Cocaine", "Sucrose"}:
        sys.exit(f"ERROR: unexpected treatment '{treatment}' in '{sample_name}'.")
    return sex, treatment, replicate


def main():
    sc.settings.seed = RANDOM_SEED
    sample_dirs = find_sample_dirs()
    print(f"Found {len(sample_dirs)} samples: {[p.name for p in sample_dirs]}")
    if len(sample_dirs) != 8:
        print(f"WARNING: expected 8 samples, found {len(sample_dirs)}. Continuing anyway.")

    adatas = {}
    for path in sample_dirs:
        name = path.name
        sex, treatment, replicate = parse_metadata(name)

        # var_names='gene_symbols' gives readable names (repo, elav, ple)
        # instead of FBgn IDs. make_unique handles the handful of symbols
        # that map to more than one FlyBase gene.
        a = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
        a.var_names_make_unique()

        a.obs["sample"] = name
        a.obs["sex"] = sex
        a.obs["treatment"] = treatment
        a.obs["replicate"] = replicate
        a.obs["condition"] = f"{sex}_{treatment}"

        print(f"  {name}: {a.n_obs:,} cells x {a.n_vars:,} genes")
        adatas[name] = a

    # join='outer' keeps every gene seen in any sample, filling absent genes
    # with zeros. join='inner' would silently drop genes missing from a single
    # sample -- including, potentially, a marker you need later.
    # index_unique='-' appends the sample key to barcodes so identical
    # barcodes from different runs do not collide.
    adata = ad.concat(adatas, label="sample_batch", index_unique="-", join="outer")
    adata.var_names_make_unique()

    # Make the categoricals ordered sensibly so plots come out predictably.
    for col in ["sample", "sex", "treatment", "replicate", "condition"]:
        adata.obs[col] = adata.obs[col].astype("category")

    print(f"\nMerged dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print("(Paper reports 86,224 cells AFTER their QC -- expect more than that here.)")
    print("\nCells per sample:")
    print(adata.obs["sample"].value_counts().sort_index().to_string())

    adata.write(H5AD_RAW)
    print(f"\nWrote {H5AD_RAW}")


if __name__ == "__main__":
    main()
