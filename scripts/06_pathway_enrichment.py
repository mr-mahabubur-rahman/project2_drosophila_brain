"""
06_pathway_enrichment.py — pathway enrichment of cocaine-responsive genes.

WHAT THIS STEP IS FOR
Research Question 4. A list of 400 gene names is not an answer; "cocaine
remodels glutathione metabolism and Toll-like receptor signalling in glia" is.
Enrichment converts the gene list into biology.

A BUG IN THE PROJECT GUIDE, WORTH KNOWING ABOUT
The guide passes gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human"]
with organism="fly". Enrichr keeps SEPARATE library catalogues per organism,
and those two names belong to the Human catalogue. Requesting them in the Fly
modality returns nothing useful -- often silently. This script queries the
live list of Fly libraries and picks from what actually exists. Flagging this
in your peer review of another group would be a genuinely useful catch.

A CAVEAT TO WRITE INTO YOUR REPORT
Enrichr's background is all annotated fly genes, but our DE test could only
ever have found genes that survived QC and are expressed in brain. Using the
full genome as background inflates significance for brain-expressed terms.
Say so; do not pretend the adjusted p-values are exact.

Requires internet (Enrichr is a web API).

Run:  uv run scripts/06_pathway_enrichment.py
"""

import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from config import (
    TABLE_DIR, FIG_DIR, ENRICHR_ORGANISM, ENRICHR_LIBRARIES_PREFERRED,
    N_GENES_FOR_ENRICHMENT, LOG2FC_THRESHOLD, PADJ_THRESHOLD,
)

try:
    import gseapy as gp
except ImportError:
    sys.exit("gseapy not installed. Run: uv pip install gseapy")


def resolve_libraries():
    """Ask Enrichr what Fly libraries exist, rather than guessing names."""
    try:
        available = gp.get_library_name(organism=ENRICHR_ORGANISM)
    except Exception as e:
        print(f"Could not reach Enrichr to list libraries ({e}).")
        print("Falling back to configured names; check results carefully.")
        return ENRICHR_LIBRARIES_PREFERRED

    chosen = [lib for lib in ENRICHR_LIBRARIES_PREFERRED if lib in available]

    # Fill gaps by matching on keyword, since Enrichr renames libraries by year.
    for keyword in ["GO_Biological_Process", "KEGG", "GO_Molecular_Function"]:
        if not any(keyword in c for c in chosen):
            match = [lib for lib in available if keyword in lib]
            if match:
                chosen.append(sorted(match)[-1])  # newest version

    if not chosen:
        print(f"No usable Fly libraries found. Available: {available[:20]}")
        sys.exit(1)

    print(f"Using Enrichr Fly libraries: {chosen}")
    return chosen


def run_enrichment(genes, label, libraries):
    """Enrichr over-representation test on one gene list."""
    genes = [g for g in dict.fromkeys(genes) if isinstance(g, str)]
    if len(genes) < 5:
        print(f"  {label}: only {len(genes)} genes, skipping.")
        return None

    print(f"  {label}: {len(genes)} genes")
    try:
        enr = gp.enrichr(gene_list=genes, gene_sets=libraries,
                         organism=ENRICHR_ORGANISM, outdir=None)
    except Exception as e:
        print(f"    Enrichr call failed: {e}")
        return None

    res = enr.results
    if res is None or res.empty:
        print("    No enrichment returned.")
        return None

    res = res.sort_values("Adjusted P-value")
    res.to_csv(TABLE_DIR / f"enrichment_{label}.csv", index=False)

    sig = res[res["Adjusted P-value"] < PADJ_THRESHOLD]
    print(f"    {len(sig)} terms at adjusted p < {PADJ_THRESHOLD}")
    if len(sig):
        for _, r in sig.head(5).iterrows():
            print(f"      {r['Term'][:62]:<62} p_adj={r['Adjusted P-value']:.2e}")
    return res


def barplot(res, label, n=12):
    """Report Figure 5: enrichment bar chart."""
    import numpy as np
    d = res.nsmallest(n, "Adjusted P-value").iloc[::-1]
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 0.36 * len(d) + 1.6))
    ax.barh(range(len(d)), -np.log10(d["Adjusted P-value"].clip(lower=1e-300)),
            color="steelblue")
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([t[:60] for t in d["Term"]], fontsize=7)
    ax.axvline(-np.log10(PADJ_THRESHOLD), ls="--", lw=0.8, c="crimson")
    ax.set_xlabel("-log10 adjusted p")
    ax.set_title(label.replace("_", " "), fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"06_enrichment_{label}.png", dpi=150)
    plt.close(fig)


def top_genes(df, direction, n):
    """
    Rank by significance, split by direction.

    Splitting up- from down-regulated matters: mixing them lets an enriched
    term hide the fact that half its genes went one way and half the other,
    which is biologically meaningless.
    """
    d = df[(df["pvals_adj"] < PADJ_THRESHOLD) &
           (df["logfoldchanges"].abs() > LOG2FC_THRESHOLD)]
    d = d[d["logfoldchanges"] > 0] if direction == "up" else d[d["logfoldchanges"] < 0]
    return d.nsmallest(n, "pvals_adj")["names"].tolist()


def main():
    libraries = resolve_libraries()

    # ---- Global lists, per sex and direction ----------------------------
    print("\n=== Global cocaine-responsive genes ===")
    for sex in ["male", "female"]:
        path = TABLE_DIR / f"de_global_{sex}_all.csv"
        if not path.exists():
            print(f"{path.name} missing -- run 05_de_analysis.py first.")
            continue
        df = pd.read_csv(path)
        for direction in ["up", "down"]:
            genes = top_genes(df, direction, N_GENES_FOR_ENRICHMENT)
            label = f"{sex}_{direction}"
            res = run_enrichment(genes, label, libraries)
            if res is not None:
                barplot(res, label)

    # ---- Per-cluster lists for the top responders ------------------------
    # This is where the paper's most interesting findings live (their C22
    # surface glia: Notch, GABA degradation, NF-kB, Toll-like receptor).
    pc_path = TABLE_DIR / "de_per_cluster_significant_genes.csv"
    summary_path = TABLE_DIR / "de_per_cluster_summary.csv"
    if pc_path.exists() and summary_path.exists():
        print("\n=== Top responding clusters ===")
        genes_df = pd.read_csv(pc_path)
        summary = pd.read_csv(summary_path)
        top_clusters = (
            summary.dropna(subset=["n_sig"])
            .groupby("cluster")["n_sig"].sum()
            .nlargest(4).index.tolist()
        )
        print(f"Analysing: {top_clusters}")
        for cl in top_clusters:
            for sex in ["Male", "Female"]:
                sub = genes_df[(genes_df["cluster"].astype(str) == str(cl)) &
                               (genes_df["sex"] == sex)]
                if len(sub) < 10:
                    continue
                genes = sub.nsmallest(N_GENES_FOR_ENRICHMENT, "pvals_adj")["names"].tolist()
                safe = str(cl).replace(" ", "-").replace("/", "-")
                label = f"cluster_{safe}_{sex.lower()}"
                res = run_enrichment(genes, label, libraries)
                if res is not None:
                    barplot(res, label)
    else:
        print("\nPer-cluster DE tables not found -- run 05_de_analysis.py first.")

    print(f"\nEnrichment tables -> {TABLE_DIR}")
    print(f"Enrichment figures -> {FIG_DIR}")
    print("\nNow compare against the paper's Supplemental Table S10. Terms to look")
    print("for: inositol phosphate metabolism (Kenyon cells), GPCR and glutamate")
    print("receptor signalling (their C16), Notch/NF-kB/TLR and glutathione")
    print("metabolism (their C22 surface glia).")


if __name__ == "__main__":
    main()
