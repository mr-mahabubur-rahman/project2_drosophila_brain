# Report outline

Structure required by §9 of the project guide. Each section lists the outputs
that feed it, so you are assembling rather than starting from a blank page.

---

## 1. Title & Abstract
Goals, design (2×2×2: sex × treatment × replicate, 8 samples), key findings.
Write this **last**.

## 2. Introduction
Cocaine neurobiology in *Drosophila*; monoamine reuptake inhibition; mushroom
body anatomy and its role as the experience-dependent learning centre; why
single-cell rather than bulk (the paper's answer: cell-type-specific responses
average out in bulk); sexual dimorphism as prior expectation from Highfill
et al. 2019 and the behavioural data in Figure 1 of the paper.

## 3. Methods
Pull the numbers from `scripts/config.py` — every parameter is there with its
rationale in a comment.

- Data source: GSE152495, 8 10x samples
- QC thresholds and cells retained → `results/tables/qc_summary_before.csv`, `..._after.csv`
- Normalization: CPM (target 1e4) + log1p; 2,000 HVGs selected per-sample
- Regression of `total_counts` and `pct_counts_mt`; scaling with max_value 10
- PCA (n_pcs), neighbours (n=15), UMAP, Leiden resolution 0.8
- DE: Wilcoxon rank-sum, `reference='Sucrose'`, thresholds |log2FC| > 1 and BH-adjusted p < 0.05
- Enrichment: Enrichr Fly libraries via GSEApy
- **State explicitly** where you diverge from the paper (Seurat/SCTransform/MAST/Bonferroni) — see README §5

## 4. Results

### Figure 1 — QC across sex and treatment
`results/figures/violin_02_qc_before.png`, `violin_02_qc_after.png`
Report cells before and after filtering, per sample.

### Figure 2 — UMAP and clustering (RQ1)
`umap_03_umap_clusters.png`, `03_resolution_sweep.png`, `umap_03_umap_metadata.png`
State your cluster count at resolution 0.8, and the resolution that lands
nearest 36. Include the UMAP coloured by sample as evidence for or against a
batch effect. Compare with Figure 2 of the paper.

### Figure 3 — Cell-type annotation
`dotplot_04_canonical_markers_dotplot.png`, `matrixplot_04_...`
Table of cluster → cell type → supporting markers, from your filled-in
`annotation_worksheet_FILLED.csv`. Note which clusters you could not annotate
confidently — the paper could not annotate their C16 either, and said so.

### Figure 4 — Differential expression (RQ2, RQ3)
`05_volcano_male.png`, `05_volcano_female.png`, `05_de_burden_heatmap.png`,
`05_sex_concordance.png`
- Male and female DE gene counts, and the ratio. Paper: 691 vs 322 (~2.1×).
- Which cell types carry the largest DE burden. Compare against their Kenyon
  cells / astrocytes / surface glia — **matched by markers, not cluster number**.
- Shared vs sex-specific genes, and those changing in opposite directions.
- Positive control: `results/tables/paper_gene_check.csv`

### Figure 5 — Pathway enrichment (RQ4)
`06_enrichment_*.png`
Compare with Supplemental Table S10: inositol phosphate metabolism in Kenyon
cells; GPCR and glutamate receptor signalling in their C16; Notch, NF-κB, TLR
and glutathione metabolism in surface glia.

### Supplementary — pseudobulk validation
`05b_pseudobulk_concordance_*.png`. Direction agreement and Spearman rho.

## 5. Discussion
- Sexual dimorphism: does your male bias match theirs in direction and rough magnitude?
- Vulnerability of Kenyon cells and glia — link to the mushroom body's role in
  reward learning and the surface glia's role as blood–brain barrier.
- **Where you differ from the paper, and why.** Methodological divergence is an
  explanation; do not present it as a failure, and do not hide it.
- **Limitations, stated plainly:**
  - Pseudoreplication — cell-level Wilcoxon inflates significance; n=2 per group
  - Enrichment background is the full annotated genome, not brain-expressed genes
  - Acute exposure only; a single time point gives no temporal information
  - The paper's own caveat: expression change is not required for behavioural effect

## 6. References
Baker et al. 2021 *Genome Res* 31:1927–1937 is the anchor. Cite Scanpy (Wolf
et al. 2018), Leiden (Traag et al. 2019), UMAP (McInnes et al. 2018), GSEApy,
and FlyBase for every marker gene you assert.

## 7. Appendix — AI Usage Disclosure
`reports/AI_USAGE_DISCLOSURE.md`. Fill it in **as you go**. Reconstructing
prompts from memory on the last day is how that section ends up inaccurate,
and inaccuracy there is an integrity problem rather than a formatting one.

---

## Honesty note

The most common failure mode in this kind of reanalysis is quietly adjusting
parameters until the numbers match the paper, then reporting agreement. If you
get 31 clusters, report 31 and explain the methodological reasons. A report
that says "we did not reproduce their exact count, and here is why" scores
better than one claiming a match that a reviewer can see was tuned into
existence.
