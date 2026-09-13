"""
04_annotate_clusters.py — find cluster markers and assign cell-type labels.

WHAT THIS STEP IS FOR
A cluster number means nothing. "Cluster 14 responds strongly to cocaine" is
not a finding; "Kenyon cells respond strongly to cocaine" is. This step turns
numbers into biology, and it is the step that requires YOUR judgement rather
than a function call.

CRITICAL WARNING ON CLUSTER NUMBERS
The paper's C11 is Kenyon cells. Your cluster 11 is almost certainly something
else -- Leiden numbers clusters by size, and your clustering is not theirs.
Match by MARKER GENES, never by number. Writing "our C11 agrees with their C11"
without checking markers is the error most likely to be caught in peer review.

This script produces Figure 3 of your report and an annotation worksheet you
must fill in by hand.

Run:  uv run scripts/04_annotate_clusters.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from config import (
    H5AD_CLUSTERED, H5AD_ANNOTATED, FIG_DIR, TABLE_DIR,
    LEIDEN_KEY, CANONICAL_MARKERS, DE_METHOD,
)

sc.settings.verbosity = 3
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)


def available(adata, genes):
    """Keep only genes actually present -- missing ones crash the plot calls."""
    return [g for g in genes if g in adata.raw.var_names]


def main():
    adata = sc.read_h5ad(H5AD_CLUSTERED)
    n_clusters = adata.obs[LEIDEN_KEY].nunique()
    print(f"Loaded {adata.n_obs:,} cells, {n_clusters} clusters")

    # ---- 1. Cluster marker genes ----------------------------------------
    # use_raw=True is essential: adata.X here is scaled z-scores restricted to
    # 2,000 HVGs. adata.raw holds log-normalized counts for ALL genes, which is
    # what a fold change should be computed on.
    #
    # The paper used |log_e FC| > 0.5 with Bonferroni correction for markers.
    # Wilcoxon is scanpy's robust default and does not assume normality --
    # appropriate for the zero-inflated, non-Gaussian counts we have here.
    sc.tl.rank_genes_groups(
        adata, groupby=LEIDEN_KEY, method=DE_METHOD, use_raw=True,
        pts=True,  # also report the fraction of cells expressing each gene
    )

    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(TABLE_DIR / "cluster_markers_all.csv", index=False)

    # Top 5 per cluster -- this is the table you actually read while annotating.
    top5 = (
        markers.sort_values(["group", "scores"], ascending=[True, False])
        .groupby("group", observed=True)
        .head(5)
    )
    top5.to_csv(TABLE_DIR / "cluster_markers_top5.csv", index=False)

    print("\nTop 3 markers per cluster:")
    for cl in sorted(adata.obs[LEIDEN_KEY].unique(), key=lambda x: int(x)):
        genes = markers[markers["group"] == cl].nlargest(3, "scores")["names"].tolist()
        n_cells = int((adata.obs[LEIDEN_KEY] == cl).sum())
        print(f"  cluster {cl:>3} (n={n_cells:>6,}): {', '.join(genes)}")

    sc.pl.rank_genes_groups_dotplot(
        adata, n_genes=3, groupby=LEIDEN_KEY, show=False,
        save="_04_top_markers_dotplot.png",
    )

    # ---- 2. Canonical marker panel --------------------------------------
    # This is the plot you read row by row to assign identities.
    # repo high + elav low  -> glia.  elav/nSyb high -> neuron.
    # Then subdivide neurons by neurotransmitter (VAChT, Gad1, VGlut, ple...).
    flat = []
    var_groups = {}
    for cell_type, genes in CANONICAL_MARKERS.items():
        present = available(adata, genes)
        if present:
            var_groups[cell_type] = present
            flat.extend(present)

    missing = [g for gs in CANONICAL_MARKERS.values() for g in gs
               if g not in adata.raw.var_names]
    if missing:
        print(f"\nMarkers not found in the dataset (check spelling/annotation): {missing}")

    sc.pl.dotplot(
        adata, var_names=var_groups, groupby=LEIDEN_KEY, use_raw=True,
        standard_scale="var", dendrogram=True, show=False,
        save="_04_canonical_markers_dotplot.png",
    )
    sc.pl.matrixplot(
        adata, var_names=var_groups, groupby=LEIDEN_KEY, use_raw=True,
        standard_scale="var", dendrogram=True, cmap="viridis", show=False,
        save="_04_canonical_markers_matrix.png",
    )

    key_umap = available(adata, ["repo", "elav", "nSyb", "ey", "Fas2",
                                 "VAChT", "Gad1", "VGlut", "ple", "SerT", "Tdc2", "ninaE"])
    sc.pl.umap(adata, color=key_umap, ncols=4, use_raw=True, show=False,
               save="_04_umap_markers.png")

    # ---- 3. Mean marker expression table --------------------------------
    # The dotplot is for the report; this CSV is for you, because reading exact
    # numbers beats squinting at dot sizes when two clusters look similar.
    raw_df = sc.get.obs_df(adata, keys=flat, use_raw=True)
    raw_df[LEIDEN_KEY] = adata.obs[LEIDEN_KEY].values
    mean_expr = raw_df.groupby(LEIDEN_KEY, observed=True).mean().round(3)
    mean_expr.to_csv(TABLE_DIR / "marker_mean_expression_by_cluster.csv")
    print(f"\nWrote marker expression table to {TABLE_DIR/'marker_mean_expression_by_cluster.csv'}")

    # ---- 4. Annotation worksheet ----------------------------------------
    # Deliberately left blank. Filling this in IS the intellectual work of the
    # step, and you must be able to defend every label at grading.
    counts = adata.obs[LEIDEN_KEY].value_counts().sort_index()
    worksheet = pd.DataFrame({
        "cluster": counts.index,
        "n_cells": counts.values,
        "pct_of_total": (100 * counts.values / adata.n_obs).round(2),
        "top_markers": [
            ", ".join(markers[markers["group"] == cl].nlargest(5, "scores")["names"])
            for cl in counts.index
        ],
        "cell_type_ANNOTATE_ME": "",
        "paper_cluster_match": "",
        "evidence_notes": "",
    })
    ws_path = TABLE_DIR / "annotation_worksheet.csv"
    worksheet.to_csv(ws_path, index=False)

    # ---- 5. Apply annotations if the worksheet has been filled in -------
    filled = TABLE_DIR / "annotation_worksheet_FILLED.csv"
    if filled.exists():
        ann = pd.read_csv(filled, dtype={"cluster": str})
        mapping = dict(zip(ann["cluster"], ann["cell_type_ANNOTATE_ME"]))
        adata.obs["cell_type"] = (
            adata.obs[LEIDEN_KEY].astype(str).map(mapping).fillna("Unannotated")
        ).astype("category")
        sc.pl.umap(adata, color="cell_type", legend_loc="on data", legend_fontsize=5,
                   show=False, save="_04_umap_celltypes.png")
        print(f"Applied annotations from {filled.name}")
    else:
        adata.obs["cell_type"] = "Unannotated"
        print(f"\n{'='*68}\nNEXT ACTION REQUIRED")
        print(f"1. Open {ws_path}")
        print("2. Fill 'cell_type_ANNOTATE_ME' using the dotplots in results/figures/")
        print("   and Supplemental Table S4 of the paper.")
        print("3. Save it as annotation_worksheet_FILLED.csv in the same folder.")
        print("4. Re-run this script to apply the labels.")
        print(f"{'='*68}")

    adata.write(H5AD_ANNOTATED)
    print(f"\nWrote {H5AD_ANNOTATED}")
    plt.close("all")


if __name__ == "__main__":
    main()
