# Methods and Abstract — draft

---

## Abstract

Acute cocaine exposure alters gene expression across the *Drosophila* brain in a
sexually dimorphic manner (Baker et al., *Genome Research* 2021). We reanalysed
the underlying single-cell data (GEO: GSE152495; 8 samples, sex × treatment ×
two replicates) with an independent Scanpy-based pipeline to test which of the
published findings reproduce.

Quality filtering retained 86,177 cells, within 0.05% of the 86,224 reported.
Leiden clustering at the published resolution yielded 30 clusters rather than
36; a resolution sweep reached 36 at approximately 1.33, and the reported
plateau in cluster count at resolution 0.8 was not observed. Twenty-four
clusters were assigned cell types using two independent marker panels, including
populations matching the paper's Kenyon cells, surface glia, and unannotated
antennal/optic lobe cluster.

The reported male bias in transcriptional response reproduced, but its magnitude
depended on the unit of analysis: 6.43-fold when cells of each sex were pooled,
against 1.65-fold when 23 clusters were tested separately and summed (published
value 2.15-fold). Female responses exceeded male in 6 of 23 clusters, including
Kenyon cells and surface glia.

The most significant genes in each sex tracked technical covariates rather than
treatment. In males, mitochondrial and respiratory-chain transcripts dominated,
and mitochondrial content differed 6.4-fold between treatment arms (median 0.373%
versus 0.058%); oxidative phosphorylation was the top enriched pathway in three
separate male clusters and in none in females, which showed no such difference in
mitochondrial content. In females, the leading genes were yolk proteins present as
ambient RNA, detected in 27.2% of cells from one cocaine replicate against 3.2–4.9%
in controls and distributed across cell types that do not transcribe them.

The design-level conclusions of the original study reproduce robustly. The
gene- and pathway-level findings appear substantially confounded by technical
variation not reported in the original analysis.

*(~290 words. Trim the third paragraph first if you need to reach 250.)*

---

## Methods

### Data

Count matrices for GEO accession GSE152495 (Baker et al. 2021) were obtained
pre-processed, comprising eight 10x Genomics samples in a 2 × 2 × 2 design (sex:
female, male; treatment: cocaine, sucrose; two biological replicates per
combination). Matrices were produced by the original authors with CellRanger
v3.1 against *D. melanogaster* Release 6 (GCA_000001215.4).

Three corrections were applied before analysis.

**Delimiter.** The supplied `features.tsv.gz` files were space- rather than
tab-delimited, which `scanpy.read_10x_mtx` cannot parse. They were converted to
tab-delimited with field content unchanged.

**Gene nomenclature.** *nanchung* (Dmel_CG5842, symbol `nan`) is interpreted by
pandas as a missing value; it was renamed `CG5842`. The reference additionally
substitutes vertebrate ortholog names for five fly symbols, resolved by gene ID:
`VGlut1` = VGlut (CG9887), `Adcy1` = *rut* (CG9533), `Pde4` = *dnc* (CG32498),
`Trhn` = Trh (CG9122), `Tret1` = Tret1-1 (CG30035). Note that `trh` (lowercase)
in this reference is CG42865, *trachealess* — an unrelated gene — so marker
matching was performed case-sensitively.

**Sample labels.** Sex-specific markers separated samples according to the
treatment field of the folder name rather than the sex field. Sample identities
were reconstructed from the deposition order given in the authors' published
analysis code and verified against two independent marker panels
(male: `roX1`, `roX2`, MSL complex; female: `Sxl`, `Yp1`–`Yp3`, `tra`). Four sex
and four treatment labels were corrected. Original labels are retained in the
object metadata.

### Quality control

Cells with fewer than 300 or more than 2,500 detected genes were removed,
following the original study's criteria, as were cells exceeding 10%
mitochondrial reads and genes detected in fewer than 5 cells. Mitochondrial
genes were identified by the `mt:` prefix and ribosomal genes by `RpS`/`RpL`.
86,177 of 88,991 cells (96.8%) and 12,189 of 17,481 genes were retained.

### Normalisation and dimensionality reduction

Counts were normalised to 10,000 per cell and log₁ₚ-transformed. The
log-normalised matrix was retained for differential expression. Two thousand
highly variable genes were selected with per-sample batch correction
(`batch_key='sample'`, seurat flavour). The matrix was subset to these genes;
`total_counts` and `pct_counts_mt` were regressed out and values scaled to unit
variance with a maximum of 10.

Principal component analysis used the ARPACK solver. The first 30 components
(16.9% of variance) were used for a 15-nearest-neighbour graph and for UMAP.
Leiden clustering was performed at resolution 0.8, with a sweep across 0.4–1.4.
All stochastic steps used a fixed random seed (0).

### Cell-type annotation

Cluster markers were identified by Wilcoxon rank-sum test against all other
cells, computed on log-normalised values. Cell types were assigned using two
independent panels: the published top-three markers for each of the original 36
clusters, and a canonical *Drosophila* neuronal and glial panel. Assignments
required agreement between panels; clusters resolving to different identities
under the two panels were left unannotated. Twenty-four of 30 clusters were
assigned; 13,526 cells (15.7%) remained unannotated.

One cluster (n = 3,196) was excluded from differential expression analysis: it
showed mixed neuronal and glial markers and a 1.50-fold difference in median
detected genes between treatment arms — the largest of any cluster — making its
299 differentially expressed genes uninterpretable.

### Differential expression

Cocaine was compared with sucrose using the Wilcoxon rank-sum test on
log-normalised values, with sucrose as reference, so that positive log₂ fold
changes indicate higher expression in cocaine. Genes were considered
differentially expressed at |log₂FC| > 1 and Benjamini–Hochberg-adjusted
*p* < 0.05. Note that the original study used a natural-log threshold
(|ln FC| > 1, equivalent to |log₂FC| > 1.44) with Bonferroni correction and the
MAST test; gene counts are therefore not directly comparable.

Three comparisons were performed: pooled across all cells within each sex;
within each cluster and sex separately, excluding clusters with fewer than 20
cells in either arm; and overlap between the sexes.

### Pseudobulk validation

To assess pseudoreplication, counts were summed within each sample to give one
profile per biological replicate, converted to log₂ counts per million, and
compared by two-sample *t*-test with Benjamini–Hochberg correction (n = 2 per
arm). Genes with fewer than 50 total counts were excluded. Concordance with the
cell-level result was assessed by sign agreement and Spearman correlation over
the genes called significant at cell level.

### Pathway enrichment

Over-representation analysis used GSEApy against the Enrichr *Drosophila*
libraries (GO Biological Process 2018, GO Molecular Function 2018, KEGG 2019).
Up- and down-regulated genes were tested separately, up to 150 genes per list
ranked by adjusted *p*-value. Terms were considered enriched at adjusted
*p* < 0.05. Enrichr tests against all annotated fly genes, whereas the
differential expression test could only detect brain-expressed genes surviving
quality control; adjusted *p*-values are therefore optimistic for
brain-expressed terms.

### Software

Python 3.11.15 with scanpy 1.11.5, anndata 0.12.19, numpy 2.4.6, pandas 2.3.3,
scipy 1.17.1, scikit-learn 1.9.1, python-igraph 1.0.0, leidenalg 0.12.0,
umap-learn 0.5.12, gseapy 1.3.1 and matplotlib 3.11.2. The environment was
managed with `uv`; a complete lockfile and all analysis code are available at
*[repository URL]*.

### Data and code availability

Raw data: GEO accession GSE152495. Analysis code, environment lockfile, result
tables and figure-generating scripts: *[repository URL]*. Count matrices are not
redistributed; scripts to rebuild them from GEO are included.

---

## Notes

- **Insert your repository URL** in two places above.
- **The Methods paragraph on the excluded cluster** is unusual to include but
  necessary — excluding data without saying so is worse than the exclusion.
- **The abstract runs ~290 words.** If a limit applies, the third paragraph
  (pooled versus per-cluster) compresses to one sentence, though it is arguably
  your most interesting result.
- **Consider adding a Supplementary Note** listing the three data corrections
  with the diagnostic output, so a reader can verify them independently.
