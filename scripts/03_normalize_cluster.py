"""
03_normalize_cluster.py — normalize, select HVGs, reduce dimensions, cluster.

WHAT THIS STEP IS FOR
This is where Research Question 1 (reproduce the 36 clusters) is answered, and
it is the slowest step. Budget 20-60 minutes depending on your RAM.

THE ONE THING TO UNDERSTAND HERE:
adata.X gets overwritten repeatedly. After scaling it holds z-scores, which are
meaningless for differential expression (you cannot take a fold change of a
z-score). So before scaling we stash the log-normalized matrix in
`adata.raw` and in `adata.layers['lognorm']`. Every DE call downstream must
use one of those, never the scaled X. This is the single most common bug in
student scRNA-seq pipelines.

Run:  uv run scripts/03_normalize_cluster.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from config import (
    H5AD_QC, H5AD_CLUSTERED, FIG_DIR, TABLE_DIR,
    TARGET_SUM, N_TOP_GENES, HVG_BATCH_KEY,
    REGRESS_OUT, REGRESS_VARS, SCALE_MAX_VALUE,
    N_PCS, N_NEIGHBORS, LEIDEN_RESOLUTION, LEIDEN_KEY, RANDOM_SEED,
    RESOLUTION_SWEEP, TARGET_N_CLUSTERS, USE_HARMONY, HARMONY_BATCH_KEY,
    SUBSAMPLE_N_CELLS,
)

sc.settings.verbosity = 3
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)


def run_leiden(adata, resolution, key, rep="X_pca"):
    """
    Leiden with a version-tolerant call.

    Scanpy >= 1.10 defaults to flavor='leidenalg' but warns, and will switch to
    the faster igraph backend. We ask for igraph explicitly and fall back
    cleanly on older scanpy so the script runs on whatever your team installed.
    """
    try:
        sc.tl.leiden(
            adata, resolution=resolution, key_added=key,
            flavor="igraph", n_iterations=2, directed=False,
            random_state=RANDOM_SEED,
        )
    except (TypeError, ValueError):
        sc.tl.leiden(
            adata, resolution=resolution, key_added=key,
            random_state=RANDOM_SEED,
        )
    return adata.obs[key].nunique()


def main():
    adata = sc.read_h5ad(H5AD_QC)
    print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ---- 0. Optional subsampling (low-RAM escape hatch) -----------------
    # Only for machines that cannot hold the full matrix. This CHANGES your
    # cluster count and DE power, so it must be declared in Methods. Prefer
    # setting REGRESS_OUT=False first -- that saves memory without discarding
    # data.
    if SUBSAMPLE_N_CELLS is not None and adata.n_obs > SUBSAMPLE_N_CELLS:
        sc.pp.subsample(adata, n_obs=SUBSAMPLE_N_CELLS, random_state=RANDOM_SEED)
        print(f"SUBSAMPLED to {adata.n_obs:,} cells -- declare this in your Methods.")

    # ---- 1. Keep the raw counts ----------------------------------------
    # Some tools (scVI, DESeq2-style pseudobulk) need integer counts back.
    # Storing them costs memory but saves re-running steps 01-02.
    adata.layers["counts"] = adata.X.copy()

    # ---- 2. Normalize ---------------------------------------------------
    # Library-size normalization: a cell with 20,000 UMIs is not "expressing
    # more" than one with 5,000, it was just sequenced deeper. Scale every cell
    # to the same total, then log1p to tame the long right tail so that a
    # 2-fold change means the same thing at high and low expression.
    #
    # The paper used SCTransform (regularized negative binomial) instead.
    # That is a genuinely different normalization and is part of why your
    # cluster count will differ from theirs. Say so in your Methods.
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)

    # float32 halves the memory of every matrix from here on. scRNA-seq counts
    # carry nowhere near float64 precision, so nothing is lost.
    adata.X = adata.X.astype("float32")

    # adata.raw holds the full gene set, log-normalized -- this is what every
    # DE call downstream uses via use_raw=True.
    # We deliberately do NOT also keep a layers['lognorm'] copy: it would be a
    # third sparse copy of the same matrix (~800 MB at this cell count) for no
    # benefit.
    adata.raw = adata

    # ---- 3. Highly variable genes ---------------------------------------
    # Most of ~13,000 genes are uninformative noise for clustering. Keeping the
    # top 2,000 most variable ones improves signal-to-noise and speed.
    # batch_key ranks genes by how consistently variable they are ACROSS
    # samples, so a gene that is wild in one bad sample does not get picked.
    sc.pp.highly_variable_genes(
        adata, n_top_genes=N_TOP_GENES, batch_key=HVG_BATCH_KEY, flavor="seurat"
    )
    sc.pl.highly_variable_genes(adata, show=False, save="_03_hvg.png")
    print(f"Selected {int(adata.var['highly_variable'].sum())} HVGs")

    # Subset to HVGs. This is the key memory/time decision: regress_out and
    # scale on 2,000 genes take minutes; on 13,000 they take hours and can
    # exhaust RAM at ~80k cells. adata.raw still holds all genes for DE.
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"Subset to HVGs: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ---- 4. Regress out technical covariates ----------------------------
    if REGRESS_OUT:
        print(f"Regressing out {REGRESS_VARS} (slow -- this is the long wait)...")
        sc.pp.regress_out(adata, REGRESS_VARS)
    else:
        print("Skipping regress_out (REGRESS_OUT=False in config.py).")

    # Z-score each gene so PCA is not dominated by a handful of very highly
    # expressed genes. max_value=10 clips extreme outliers.
    sc.pp.scale(adata, max_value=SCALE_MAX_VALUE)

    # ---- 5. PCA ----------------------------------------------------------
    sc.tl.pca(adata, svd_solver="arpack", use_highly_variable=True, random_state=RANDOM_SEED)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False, save="_03_pca_variance.png")
    print("\nCheck figures/pca_variance_03_pca_variance.png: the elbow tells you")
    print(f"how many PCs carry signal. config.py currently uses N_PCS={N_PCS}.")

    use_rep = "X_pca"
    if USE_HARMONY:
        try:
            sc.external.pp.harmony_integrate(adata, key=HARMONY_BATCH_KEY)
            use_rep = "X_pca_harmony"
            print("Harmony integration applied.")
        except Exception as e:
            print(f"Harmony unavailable ({e}); continuing on uncorrected PCA.")

    # ---- 6. Neighbour graph and UMAP -------------------------------------
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS,
                    use_rep=use_rep, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)

    # ---- 7. Resolution sweep (answers RQ1 properly) ----------------------
    # The paper justified resolution 0.8 by showing cluster count plateaus
    # there. Reproducing that curve is a much stronger result than reporting
    # a single number, and it is the honest way to handle "we got 29, not 36".
    print("\n--- Resolution sweep ---")
    sweep = []
    for res in RESOLUTION_SWEEP:
        key = f"leiden_sweep_{res}"
        n = run_leiden(adata, res, key)
        sweep.append({"resolution": res, "n_clusters": n})
        print(f"  resolution {res}: {n} clusters")

    sweep_df = pd.DataFrame(sweep)
    sweep_df.to_csv(TABLE_DIR / "resolution_sweep.csv", index=False)

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(sweep_df["resolution"], sweep_df["n_clusters"], "o-")
    ax.axhline(TARGET_N_CLUSTERS, ls="--", c="crimson",
               label=f"paper: {TARGET_N_CLUSTERS} clusters")
    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("Number of clusters")
    ax.set_title("Cluster count vs. resolution")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "03_resolution_sweep.png", dpi=150)
    plt.close(fig)

    best = sweep_df.iloc[(sweep_df["n_clusters"] - TARGET_N_CLUSTERS).abs().argmin()]
    print(f"\nClosest to {TARGET_N_CLUSTERS} clusters: resolution "
          f"{best['resolution']} -> {int(best['n_clusters'])} clusters")

    # ---- 8. The primary clustering (paper's resolution) ------------------
    n_main = run_leiden(adata, LEIDEN_RESOLUTION, LEIDEN_KEY)
    print(f"\nPrimary clustering at resolution {LEIDEN_RESOLUTION}: "
          f"{n_main} clusters (paper target: {TARGET_N_CLUSTERS})")

    # ---- 9. Diagnostic UMAPs ---------------------------------------------
    # The paper reports no cluster was dominated by one sample/sex/replicate.
    # These plots are how you check whether that holds for YOUR embedding.
    # If clusters segregate by sample, you have a batch effect -> set
    # USE_HARMONY = True and re-run.
    sc.pl.umap(adata, color=[LEIDEN_KEY], legend_loc="on data", legend_fontsize=6,
               show=False, save="_03_umap_clusters.png")
    sc.pl.umap(adata, color=["sex", "treatment", "condition", "sample"], ncols=2,
               show=False, save="_03_umap_metadata.png")
    sc.pl.umap(adata, color=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
               ncols=3, show=False, save="_03_umap_qc.png")

    # Quantify the "no cluster is dominated" claim instead of eyeballing it.
    comp = (
        pd.crosstab(adata.obs[LEIDEN_KEY], adata.obs["sample"], normalize="index")
        .round(3)
    )
    comp.to_csv(TABLE_DIR / "cluster_composition_by_sample.csv")
    max_frac = comp.max(axis=1)
    print(f"\nMost sample-skewed cluster: {max_frac.idxmax()} "
          f"({max_frac.max():.1%} from a single sample; even mixing would be 12.5%)")

    adata.write(H5AD_CLUSTERED)
    print(f"\nWrote {H5AD_CLUSTERED}")
    plt.close("all")


if __name__ == "__main__":
    main()
