"""
02_qc_filter.py — compute QC metrics, plot them, then filter.

WHAT THIS STEP IS FOR
Droplet data contains three things that are not cells: empty droplets with
ambient RNA (few genes), dying cells (high mitochondrial fraction), and
doublets (two cells in one droplet, so too many genes). Each distorts
clustering in a different direction.

IMPORTANT: this script PLOTS BEFORE IT CUTS. Look at
results/figures/violin_02_qc_before.png before you accept the thresholds in
config.py. Numbers copied from a tutorial are not a justification you can
defend at grading.

Produces Figure 1 of your report (QC metrics across sex and treatment).

Run:  uv run scripts/02_qc_filter.py
"""

import matplotlib
matplotlib.use("Agg")  # write files, don't try to open windows (WSL has no display)
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from config import (
    H5AD_RAW, H5AD_QC, FIG_DIR, TABLE_DIR,
    MIN_GENES_PER_CELL, MAX_GENES_PER_CELL, MIN_CELLS_PER_GENE, MAX_PCT_MT,
)

sc.settings.verbosity = 3
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)


def annotate_qc_genes(adata):
    """
    Flag mitochondrial and ribosomal genes.

    Drosophila naming is NOT the human convention -- no 'MT-' prefix.
    Mito genes are 'mt:' (mt:CoI, mt:ND1, mt:Cytb) and ribosomal proteins
    are 'RpS'/'RpL'. Getting this wrong gives every cell 0% mt and you would
    filter nothing, which is the most common silent failure in this dataset.
    """
    adata.var["mt"] = adata.var_names.str.startswith("mt:")
    adata.var["ribo"] = adata.var_names.str.startswith(("RpS", "RpL", "rps", "rpl"))

    n_mt, n_ribo = int(adata.var["mt"].sum()), int(adata.var["ribo"].sum())
    print(f"Flagged {n_mt} mitochondrial and {n_ribo} ribosomal genes.")
    if n_mt == 0:
        print(
            "WARNING: zero mitochondrial genes flagged. Your var_names are probably\n"
            "         FlyBase IDs (FBgn...) rather than symbols. Re-run step 01."
        )
    return adata


def qc_summary_table(adata, label):
    """Per-sample QC summary -- goes straight into the Methods section."""
    g = adata.obs.groupby("sample", observed=True)
    tbl = pd.DataFrame({
        "n_cells": g.size(),
        "median_genes": g["n_genes_by_counts"].median(),
        "median_counts": g["total_counts"].median(),
        "median_pct_mt": g["pct_counts_mt"].median(),
        "median_pct_ribo": g["pct_counts_ribo"].median(),
    }).round(2)
    path = TABLE_DIR / f"qc_summary_{label}.csv"
    tbl.to_csv(path)
    print(f"\nQC summary ({label}):\n{tbl.to_string()}")
    print(f"Wrote {path}")
    return tbl


def main():
    adata = sc.read_h5ad(H5AD_RAW)
    print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    adata = annotate_qc_genes(adata)

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )

    # ---- LOOK BEFORE YOU CUT -------------------------------------------
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
        groupby="condition", jitter=0.4, multi_panel=True, rotation=45,
        show=False, save="_02_qc_before.png",
    )
    sc.pl.scatter(
        adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt",
        show=False, save="_02_counts_vs_genes.png",
    )
    qc_summary_table(adata, "before")

    # ---- FILTER --------------------------------------------------------
    n0 = adata.n_obs
    print("\n--- Filtering ---")

    sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
    print(f"After min_genes >= {MIN_GENES_PER_CELL}: {adata.n_obs:,} cells "
          f"({n0 - adata.n_obs:,} removed)")

    if MAX_GENES_PER_CELL is not None:
        # Upper bound = crude doublet filter. The paper uses 2500 on Seurat's
        # advice: droplets with more detected genes than any single fly neuron
        # plausibly expresses are usually two cells stuck together.
        n_before = adata.n_obs
        adata = adata[adata.obs["n_genes_by_counts"] < MAX_GENES_PER_CELL].copy()
        print(f"After n_genes < {MAX_GENES_PER_CELL}: {adata.n_obs:,} cells "
              f"({n_before - adata.n_obs:,} removed as likely doublets)")

    n_before = adata.n_obs
    adata = adata[adata.obs["pct_counts_mt"] < MAX_PCT_MT].copy()
    print(f"After pct_mt < {MAX_PCT_MT}: {adata.n_obs:,} cells "
          f"({n_before - adata.n_obs:,} removed as stressed/dying)")

    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    print(f"After min_cells >= {MIN_CELLS_PER_GENE}: {adata.n_vars:,} genes "
          f"({n_genes_before - adata.n_vars:,} removed)")

    print(f"\nRetained {adata.n_obs:,} / {n0:,} cells ({100*adata.n_obs/n0:.1f}%)")
    print("Paper analysed 86,224 cells. If you are far off, revisit the thresholds.")

    # ---- Confirm the design survived filtering -------------------------
    # If one arm lost most of its cells, the DE comparison is compromised and
    # you must say so in the report rather than pretend the design is balanced.
    print("\nCells per condition after QC:")
    print(adata.obs["condition"].value_counts().sort_index().to_string())

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        groupby="condition", jitter=0.4, multi_panel=True, rotation=45,
        show=False, save="_02_qc_after.png",
    )
    qc_summary_table(adata, "after")

    adata.write(H5AD_QC)
    print(f"\nWrote {H5AD_QC}")
    plt.close("all")


if __name__ == "__main__":
    main()
