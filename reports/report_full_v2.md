# Independent reanalysis of the cocaine-exposed *Drosophila* brain: the atlas reproduces, but sex markers do not agree with the deposited sample labels

**Mahabubur Rahman**
Student Project 2 — scRNA-seq reanalysis of GSE152495
September 2026
Repository: https://github.com/mr-mahabubur-rahman/project2_drosophila_brain

---

## Abstract

Acute cocaine consumption is reported to alter gene expression across the
*Drosophila melanogaster* brain in a cell-type-specific and sexually dimorphic
manner (Baker et al. 2021). We reanalysed the underlying single-cell data
(GEO: GSE152495; eight samples, sex × treatment × two biological replicates)
using an independent Scanpy-based pipeline.

Quality filtering retained **86,177 cells against the 86,224 reported**, a
difference of 0.05%. Leiden clustering at the published resolution yielded 30
clusters rather than 36; a resolution sweep reached 36 at approximately 1.33, and
the reported plateau in cluster number at resolution 0.8 was not observed.
Twenty-four clusters were assigned cell-type identities from two independent
marker panels, recovering populations corresponding to Kenyon cells, surface
glia and the paper's unannotated antennal/optic lobe cluster.

Three errors in the supplied data were identified and corrected before analysis:
space-delimited `features.tsv.gz` files that prevent gene-symbol assignment; the
gene *nanchung* (`nan`), which pandas parses as a missing value; and five fly
symbols replaced by vertebrate ortholog names in the reference, one of which
(`trh` = *trachealess*) would cause a cluster to be misidentified as serotonergic
under case-insensitive matching.

**Differential expression could not be interpreted as a treatment effect.**
Sex-specific markers do not agree with the deposited sex labels: `roX1` and
`roX2` read 3.85–3.97 and 2.16–2.40 in all four sucrose samples and 0.06–0.12 and
0.008–0.015 in all four cocaine samples, while `Yp1`–`Yp3` show the reverse. The
split follows treatment, crossing both declared sexes. Consequently, when cocaine
is compared with sucrose within each declared sex, the sex markers appear among
the most significant genes in **both** comparisons — `roX1` and `roX2` at
approximately −9.5 log₂FC with adjusted *p* below machine precision.

Sample identities were verified three ways: per-sample cell counts in
Supplemental Table S2 match all eight deposited matrices exactly (total 88,991);
the GEO sample titles agree; and the authors' Supplemental Code documents the
same assignment from CellRanger output directory through to the Seurat objects.
The discrepancy therefore did not arise during analysis or deposition.

We report the reproduction results, document the discrepancy in full, and refrain
from interpreting any differential expression result as a cocaine effect. The
question has been referred to the original authors.

![](results/figures/fig0_graphical_abstract.png)

**Figure 1. Graphical abstract.** Overview of the data, pipeline, what
reproduces, and the sex-marker discrepancy.

---

## 1. Introduction

### 1.1 *Drosophila* as a model for psychostimulant response

Cocaine acts on the vertebrate nervous system principally by blocking monoamine
reuptake, elevating synaptic dopamine, serotonin and noradrenaline. The
molecular targets of this action are conserved well beyond vertebrates.
*Drosophila melanogaster* possesses a dopamine transporter whose crystal
structure accommodates cocaine (Wang et al. 2015), a cocaine-sensitive serotonin
transporter (Corey et al. 1994), and a vesicular monoamine transporter whose
overexpression attenuates behavioural responses to the drug (Chang et al. 2006).
Roughly 75% of human disease-associated genes have fly orthologs (Pandey and
Nichols 2011).

Behaviourally, flies exposed to cocaine show impaired negative geotaxis, reduced
startle responses, increased grooming and, at higher exposures, seizures
(McClung and Hirsh 1998; Filošević et al. 2018). Repeated intermittent exposure
produces sensitisation. These phenotypes are sexually dimorphic: male flies show
greater locomotor impairment and more pronounced compulsive grooming than
females (Baker et al. 2021), paralleling sex differences in human psychostimulant
use, where women escalate to compulsive use more rapidly and metabolise cocaine
faster than men (Becker and Koob 2016; Lukas et al. 1996).

The genetic architecture of cocaine consumption in flies has been mapped using
the *Drosophila* Genetic Reference Panel, implicating dopaminergic neurons and
the mushroom bodies — the central integrative structures underlying
experience-dependent modification of behaviour (Highfill et al. 2019). Notably,
RNAi knockdown in glia also altered consumption, indicating that the response is
not confined to neurons.

### 1.2 Why single-cell resolution

Bulk transcriptomic profiling of whole brain averages expression across every
cell type present. A response confined to one population — Kenyon cells, say, or
surface glia — is diluted by the tens of thousands of cells that do not respond,
and may fall below detection entirely. Worse, opposing responses in different
cell types can cancel, producing an apparent absence of effect where substantial
regulation is occurring.

Single-cell RNA sequencing removes this averaging. The adult fly brain is small
enough (~100,000 cells) that the entire organ can be surveyed in a single
experiment, and reference atlases of the aging and midbrain transcriptome
already exist (Davie et al. 2018; Croset et al. 2018). This makes *Drosophila*
unusually well suited to asking which cell populations respond to an acute drug
exposure, rather than merely whether the brain as a whole does.

### 1.3 The dataset reanalysed here

Baker et al. (2021, *Genome Research* 31:1927–1937) applied this approach to
acute cocaine consumption. Flies of both sexes consumed a fixed volume of
sucrose or cocaine-supplemented sucrose in a capillary feeder assay; brains were
dissected immediately and dissociated for 10x Genomics sequencing, in duplicate,
giving eight samples in a 2 × 2 × 2 design. Unsupervised clustering of 86,224
cells yielded 36 populations spanning the major neuronal and glial classes.

Their principal conclusions were that acute cocaine elicits rapid, widespread
transcriptional change across the brain; that the response is cell-type
specific, with Kenyon cells and glia among the strongest responders; and that it
is profoundly sexually dimorphic, with 691 differentially expressed genes in
males against 322 in females and little overlap between the sexes.

### 1.4 Rationale for independent reanalysis

Deposited single-cell data can be reanalysed independently, and doing so serves
two distinct purposes that are worth separating.

The first is **reproduction**: does an independent pipeline, using different
software and different statistical choices, recover the published result? A
finding that survives Scanpy as well as Seurat, Leiden as well as shared
nearest-neighbour clustering, and Wilcoxon as well as MAST is more robust than
one demonstrated only once.

The second is **verification**: are the published conclusions supported by the
data as analysed, independent of whether they replicate? These come apart more
often than is comfortable. Single-cell differential expression is vulnerable to
several confounds that are not always reported — pseudoreplication, where
thousands of cells from two animals are treated as thousands of independent
observations (Squair et al. 2021); ambient RNA from lysed cells partitioning
into droplets (Young and Behjati 2020); and systematic differences in capture
efficiency or cell viability between experimental groups, which produce
expression differences that are technical in origin but statistically
indistinguishable from biological ones.

None of these is detectable from a published figure. All are detectable from the
deposited data, and each requires a specific check.

### 1.5 Aims

This reanalysis therefore set out to:

1. reconstruct the processing pipeline independently in Scanpy and determine
   which published findings reproduce;
2. assign cell-type identities using marker evidence rather than cluster
   correspondence, so that comparison with the original is made on biological
   rather than numerical grounds;
3. test the differential expression results against explicit controls —
   replicate-level pseudobulk aggregation, a formal sex × treatment interaction
   model, cluster-size matching, and comparison against alternative groupings of
   the same samples; and
4. report clearly which conclusions are supported by these data, which are
   confounded, and which the design cannot resolve.

---

## 2. Methods

![](results/figures/fig0b_flowchart.png)

**Figure 2. Analysis workflow.** Numbers in brackets refer to the pipeline
scripts in `scripts/` and `tools/`.

### 2.1 Data and sample identity

Count matrices for GSE152495 were obtained pre-processed as eight CellRanger v3
output directories, produced by the original authors against *D. melanogaster*
Release 6. Three problems were corrected before analysis; each would have
propagated silently.

**Delimiter.** The supplied `features.tsv.gz` files were space- rather than
tab-delimited, which `scanpy.read_10x_mtx` cannot parse. All three fields
collapsed into one column and gene-symbol assignment failed. Files were converted
with field content unchanged (17,481 lines each); originals retained.

**Gene nomenclature.** *nanchung* (Dmel_CG5842, FlyBase symbol `nan`) is
interpreted by pandas as a missing value rather than text, leaving one gene
unnamed and preventing HDF5 serialisation. It was renamed `CG5842`, which
survives CSV round-trips. *narrow abdomen* (`na`, Dmel_CG1517) is unaffected, as
lowercase `na` is not in the default missing-value list.

The reference substitutes vertebrate ortholog names for five fly symbols,
resolved by gene ID against FlyBase (Thurmond et al. 2019): `VGlut1` = VGlut
(CG9887), `Adcy1` = *rut* (CG9533), `Pde4` = *dnc* (CG32498), `Trhn` = Trh
(CG9122), `Tret1` = Tret1-1 (CG30035). Note that `trh` (lowercase) in this
reference is Dmel_CG42865, *trachealess* — an unrelated tracheal transcription
factor, not tryptophan hydroxylase. Marker matching was therefore performed
case-sensitively; a case-insensitive match returns *trachealess* and would label
a cluster serotonergic on the strength of a tracheal gene.

**Sample identity.** Sample assignment was verified against three independent
sources, which agree:

1. **Supplemental Table S2** lists per-sample cell counts. All eight match the
   deposited matrices exactly (9,072 / 11,693 / 13,193 / 11,033 / 13,072 / 9,367
   / 10,437 / 11,124; total 88,991).
2. **GEO sample titles** correspond to the same assignment.
3. **The authors' Supplemental Code** documents the mapping from CellRanger
   output directories to Seurat objects: S1 and S2 female sucrose, S3 and S4 male
   sucrose, S5 and S6 female cocaine, S7 and S8v2 male cocaine.

Labels were therefore used as deposited. An earlier version of this analysis
applied a marker-based reassignment; Supplemental Table S2 showed that
reassignment to be incorrect, and it was reverted. The reversion is recorded in
the repository history and in `docs/data_provenance.md`.

The sex markers nonetheless disagree with these labels (§3.3), a discrepancy
that could not be resolved from the deposited data.

### 2.2 Quality control

Cells with fewer than 300 or more than 2,500 detected genes were removed,
following the original study's criteria, as were cells exceeding 10%
mitochondrial reads and genes detected in fewer than 5 cells. The original study
applied no mitochondrial threshold; it was retained here because high
mitochondrial fraction marks damaged cells, and brain dissociation is harsh.
Mitochondrial genes were identified by the `mt:` prefix and ribosomal genes by
`RpS`/`RpL`.

### 2.3 Normalisation, dimensionality reduction and clustering

Counts were normalised to 10,000 per cell and log₁ₚ-transformed. The
log-normalised matrix was retained for differential expression. Two thousand
highly variable genes were selected per sample (`batch_key='sample'`, Seurat
flavour) and ranked by consistency across samples. The matrix was subset to
these genes; `total_counts` and `pct_counts_mt` were regressed out and values
scaled to unit variance with a maximum of 10.

Principal component analysis used the ARPACK solver. The first 30 components
(16.9% of variance; Supplementary Figure S2) were used to build a
15-nearest-neighbour graph and for UMAP embedding (McInnes et al. 2018). Leiden
clustering (Traag et al. 2019) was performed at resolution 0.8, with a sweep
across 0.4–1.4. No batch integration was applied. All stochastic steps used a
fixed random seed (0).

### 2.4 Cell-type annotation

Cluster markers were identified by Wilcoxon rank-sum test of each cluster against
all other cells, on log-normalised values. Cell types were assigned using two
independent panels: the published top-three markers for each of the original 36
clusters, and a canonical *Drosophila* neuronal and glial panel with every symbol
verified by gene ID. Assignments required agreement between panels; clusters
resolving to different identities under the two panels were left unannotated.

One cluster (3,196 cells) was excluded from differential expression analysis
because it showed mixed neuronal and glial markers, so its differentially
expressed genes could not be assigned to any cell type. The exclusion was first
made under the superseded label assignment, when the cluster also showed the
largest difference in sequencing depth between treatment arms. Under the
deposited labels that difference is modest (median 745 versus 608 detected
genes per cell, 1.23-fold, against up to 1.47-fold in retained clusters), so the
exclusion now rests on mixed identity alone. Including the cluster would raise
the summed per-cluster totals in §3.4 from 1,282 and 1,256 to 1,378 and 1,308
(ratio 1.05 rather than 1.02), which does not affect any conclusion.
Differential expression was therefore carried out on 23 cell-type groups: 19
annotated cell types and four unannotated groups.

### 2.5 Differential expression

Cocaine was compared with sucrose within each declared sex using the Wilcoxon
rank-sum test on log-normalised values, with sucrose as reference, so that
positive log₂ fold changes indicate higher expression in the cocaine-labelled
samples. Genes were considered differentially expressed at |log₂FC| > 1 and
Benjamini–Hochberg-adjusted *p* < 0.05. The original study used MAST with a
natural-log threshold (|ln FC| > 1, equivalent to |log₂FC| > 1.443) and
Bonferroni correction; counts at that threshold are also reported, but gene
counts are not directly comparable between the two analyses.

Three comparisons were performed: pooled across all cells within each declared
sex; within each cell-type group and declared sex, excluding groups with fewer
than 20 cells in either arm; and overlap between the two declared sexes.

### 2.6 Sex-marker check

Mean log-normalised expression of `roX1`, `roX2` (male-specific lncRNAs of the
dosage-compensation complex), `Yp1`–`Yp3` (female yolk proteins) and `Sxl` was
computed per sample after quality control. `msl-2`, `mle` and `tra`, which are
sex-specifically regulated after transcription and therefore present in both
sexes at the RNA level, served as negative controls. The check is reproduced by
`tools/make_sexmarker_figure.py`.

### 2.7 Replicate-level controls

**Pseudobulk aggregation.** Raw counts for the 2,000 highly variable genes were
summed within each sample, giving one profile per biological replicate; genes
with fewer than 50 total counts were removed (1,755 retained). Profiles were
converted to log₂ counts per million and compared within each declared sex by
two-sample *t*-test with Benjamini–Hochberg correction (*n* = 2 versus 2).
Concordance with the cell-level result was assessed by sign agreement and
Spearman correlation over the cell-level significant genes present in this
gene set (37 in each arm).

**Interaction model.** Counts for all 12,189 post-QC genes were summed per
sample, and genes with at least 50 total counts (9,397) were converted to log₂
CPM. For each gene, an ordinary least-squares model
`expression ~ sex + treatment + sex:treatment` was fitted across the eight
samples using effect coding (±0.5), so that main effects are averages across the
other factor. Interaction and main-effect *p*-values were
Benjamini–Hochberg-adjusted separately.

### 2.8 Permutation control

Within each declared sex, the four samples can be split two-versus-two in three
distinct ways: the treatment split, the replicate axis (Cocaine R1 + Sucrose R1
versus the rest) and a diagonal pairing (Cocaine R1 + Sucrose R2 versus the
rest). The two non-treatment splits carry no treatment difference and serve as
nulls. Each was tested with the same cells, test and thresholds as the treatment
contrast.

### 2.9 Cluster-size matching

Because cluster size correlated with differentially expressed gene count
(Spearman correlation across the 23 cell-type groups, with cell and gene counts
summed across both declared sexes), each group was downsampled to 400 cells (200
per treatment arm, both declared sexes pooled) and re-tested over three random
seeds. Groups with fewer than 50 cells in either arm were skipped.

### 2.10 Pathway enrichment

Over-representation analysis used GSEApy (Fang et al. 2023) against the Enrichr
*Drosophila* libraries (Kuleshov et al. 2016): GO Biological Process 2018, GO
Molecular Function 2018 and KEGG 2019. Up- and down-regulated genes were tested
separately, with up to 150 genes per list ranked by adjusted *p*-value. Terms
were considered enriched at adjusted *p* < 0.05. Enrichr tests against all
annotated fly genes, whereas the differential expression test could only detect
brain-expressed genes surviving quality control; adjusted *p*-values are
therefore optimistic for brain-expressed terms.

### 2.11 Software and availability

Python 3.11 with scanpy 1.11.5 (Wolf et al. 2018), anndata 0.12.19, numpy 2.4.6,
pandas 2.3.3, scipy 1.17.1, scikit-learn 1.9.1, python-igraph 1.0.0, leidenalg
0.12.0, umap-learn 0.5.12, gseapy 1.3.1 and matplotlib 3.11.2. The environment
was managed with `uv`; the lockfile (`requirements.lock.txt`), all analysis code,
result tables and figure-generating scripts are available at
https://github.com/mr-mahabubur-rahman/project2_drosophila_brain. Count matrices
are not redistributed; `tools/download_data.sh` and `tools/organize_data.sh`
rebuild them from GEO.

---

## 3. Results

### 3.1 Quality control reproduces the published cell count

After filtering, **86,177 of 88,991 cells were retained (96.8%)**, against the
86,224 reported — a difference of 47 cells, or 0.05%. Gene count fell from
17,481 to 12,189. The agreement indicates the filtering criteria and input data
were correctly identified (Figure 3).

Sequencing depth varied approximately two-fold between samples (median 675–1,195
genes per cell), and was not confounded with treatment: within each declared sex,
one deeper and one shallower replicate fell in each arm.

**Mitochondrial fraction varies by sample, not by group.** Per-sample medians
span 19-fold and scatter across all four groups:

| Sample | Median mt % |
|---|---|
| Female_Cocaine_R2 | 0.700 |
| Male_Sucrose_R2 | 0.472 |
| Male_Sucrose_R1 | 0.319 |
| Male_Cocaine_R1 | 0.277 |
| Male_Cocaine_R2 | 0.251 |
| Female_Cocaine_R1 | 0.158 |
| Female_Sucrose_R1 | 0.086 |
| Female_Sucrose_R2 | 0.037 |

Cell-level group means suggest a female-specific treatment difference (0.571%
cocaine versus 0.103% sucrose) and a male difference in the opposite direction
(0.348% versus 0.528%), but both are driven by single replicates —
`Female_Cocaine_R2` and `Male_Sucrose_R2` respectively. This is
dissection-to-dissection variation, not a group-level property, and any
mitochondrial gene appearing in a treatment comparison should be interpreted with
that in mind.

![](results/figures/fig1_qc_before.png)

![](results/figures/fig1_qc_after.png)

**Figure 3. Quality control.** Genes per cell, UMI counts per cell and
mitochondrial fraction per sample (A) before filtering (88,991 cells) and (B)
after filtering (86,177 cells). Dashed lines mark the thresholds; samples are
coloured by treatment.

### 3.2 Clustering recovers 30 populations; 24 were assignable

Leiden clustering at resolution 0.8 yielded 30 clusters against the 36 reported
(Figure 4A). A resolution sweep gave a monotonic increase (19, 24, 30, 33, 34, 37
at resolutions 0.4 to 1.4), reaching 36 at approximately 1.33 (Figure 4B).
**The published stability argument — that cluster count plateaus at resolution
0.8 — was not reproduced.** This may be a property of Leiden rather than a
disagreement about the data.

The embedding is not driven by batch (Figure 4C, D). Only one of 30 clusters
exceeded 30% composition from a single sample (34.3%, against 12.5% under even
mixing), and that cluster remained 74.5% female among depth-matched cells against
a 56.4% background, indicating a genuine population rather than an artefact.

Twenty-four clusters received cell-type labels using two independent marker
panels (Figure 5). Six clusters (13,526 cells, 15.7%) were left unannotated where
the panels disagreed. Populations corresponding to the paper's Kenyon cells
(C11), surface glia (C22) and unannotated antennal/optic lobe cluster (C16) were
identified by marker gene rather than cluster number. The paper's C17
(astrocytes) could not be confidently identified: the cluster matching on `Eaat1`
was `repo`- and `alrm`-negative, while the genuine astrocyte-like population
contains only 143 cells.

![](results/figures/fig2a_umap_clusters.png)

![](results/figures/fig2b_resolution_sweep.png)

![](results/figures/fig2e_umap_by_sample.png)

![](results/figures/fig2f_cluster_composition.png)

**Figure 4. Clustering.** (A) UMAP of 86,177 cells coloured by Leiden cluster at
resolution 0.8 (30 clusters). (B) Cluster number across resolutions 0.4–1.4; the
dashed line marks the 36 clusters reported by Baker et al. (C) UMAP split by
sample. (D) Sample composition of each cluster; the dashed line marks even mixing
(12.5%).

![](results/figures/fig3a_umap_celltypes.png)

![](results/figures/fig3b_marker_dotplot.png)

**Figure 5. Cell-type annotation.** (A) UMAP coloured by assigned cell type.
(B) Expression of canonical marker genes across the 30 clusters; dot size shows
the fraction of cells expressing each gene and colour the scaled mean expression.

### 3.3 Sex markers do not agree with the deposited sample labels

This is the principal finding and it precludes interpretation of the
differential expression results below.

Five transcript-level sex markers separate the eight samples cleanly, but the
split follows the **treatment** field of the sample name rather than the **sex**
field (Figure 6; mean log-normalised expression per sample):

| Marker | Specificity | All four Sucrose samples | All four Cocaine samples |
|---|---|---|---|
| `roX1` | male (dosage-compensation lncRNA) | 3.846 – 3.972 | 0.059 – 0.116 |
| `roX2` | male (dosage-compensation lncRNA) | 2.160 – 2.396 | 0.008 – 0.015 |
| `Yp1` | female (yolk protein) | 0.000 – 0.001 | 0.065 – 0.255 |
| `Yp2` | female (yolk protein) | 0.003 – 0.003 | 0.028 – 0.153 |
| `Yp3` | female (yolk protein) | 0.006 – 0.010 | 0.069 – 0.428 |
| `Sxl` | female (supporting only) | 0.695 – 0.984 | 1.227 – 1.762 |

`roX1` and `roX2` are lncRNAs of the dosage-compensation complex, expressed
almost exclusively in males. Yolk proteins are transcribed in female fat body and
are absent from males. The ranges do not overlap for any of these five markers,
and the highest yolk-protein value in the dataset (`Yp3` at 0.428) occurs in a
sample labelled male. `Sxl` follows the same split, but because its sex
specificity arises largely through alternative splicing, which 3′ counting cannot
resolve, it is treated as supporting rather than primary evidence.

Genes whose expression is post-transcriptionally regulated and therefore
present in both sexes at the RNA level — `msl-2` (0.039–0.080), `mle`
(0.252–0.376), `tra` (0.045–0.102) — show no such separation, as expected. Only
the markers that are genuinely sex-specific at the transcript level split the
samples, which argues against a generic technical artefact.

**The consequence is direct.** Comparing cocaine with sucrose *within each
declared sex* returns the sex markers among the most significant genes in both
comparisons:

| Gene | Male-labelled contrast (log₂FC) | Female-labelled contrast (log₂FC) |
|---|---|---|
| `roX2` | −9.620 | −9.629 |
| `roX1` | −9.539 | −9.125 |
| `Yp1` | +7.684 | +6.453 |
| `Yp3` | +5.234 | +3.479 |
| `Yp2` | +5.217 | +3.829 |
| `Sxl` | +1.065 | +1.588 |

All six appear among the significant genes in both arms, with adjusted *p* at or
below machine precision. Each nominal treatment contrast is separating males from
females. The original report is consistent with this: Baker et al. list `roX2`
among the genes that responded globally to cocaine, down-regulated after
consumption (their Supplemental Tables S5 and S6).

Sample identity was verified against Supplemental Table S2, the GEO titles and
the authors' Supplemental Code (§2.1); all three agree, and no computational step
in either their pipeline or ours could have introduced the discrepancy. Neither
their code nor ours, as originally written, checked sex markers against the
sample labels.

Two explanations are consistent with the evidence and cannot be distinguished
from the deposited data: the sex annotation may have been transposed at some
point before sequencing, or sex and treatment may covary for a reason not
described in the methods. The question has been referred to the original authors.

![](results/figures/fig_sexmarker_discrepancy.png)

**Figure 6. Sex-specific markers do not agree with the deposited sample labels.**
(A) Mean log-normalised expression of male (`roX1`, `roX2`) and female (`Yp1`–`Yp3`)
markers in each sample. (B) Samples plotted by `roX1` against `Yp1` expression
(log scales) separate by declared treatment, not by declared sex. (C) Log₂ fold changes of
the sex markers in the nominal cocaine-versus-sucrose contrast within each
declared sex.

### 3.4 Differential expression results, reported but not interpreted

The results below are recorded for completeness and as evidence bearing on §3.3.
**None should be read as a cocaine effect**, since each contrast is confounded
with sex.

Pooled within each declared sex, the male-labelled contrast returned **90
differentially expressed genes** (75 up, 15 down) and the female-labelled
contrast **152** (126 up, 26 down) at |log₂FC| > 1 and BH-adjusted *p* < 0.05
(Figure 7A, B) — a ratio of 0.59, favouring the female-labelled arm. At the
paper's stricter threshold (|log₂FC| > 1.443) the counts were 20 and 61. Baker et
al. report 691 and 322, a ratio of 2.15 favouring males. Summing across the 23
cell-type groups analysed (§2.4) gives 1,282 and 1,256, a ratio of 1.02
(Figure 7C).

Overlap between the two arms was 53 genes unique to the male-labelled contrast,
115 unique to the female-labelled contrast, and 37 shared, of which one changed
direction. The shared set includes `roX1`, `roX2`, `Yp1`–`Yp3` and `Sxl`.

Cluster size correlated with differentially expressed gene count (Spearman
ρ = 0.634 across the 23 groups), so the raw ranking partly measures statistical
power. Downsampling every group to 400 cells over three seeds gave central brain
B cholinergic neurons highest (69.0 ± 8.0), followed by the unannotated
Eaat1-high cluster (43.7 ± 11.3), the antennal/optic lobe cluster (39.7 ± 15.7)
and surface glia (38.3 ± 7.4); Kenyon cells ranked sixth (26.3 ± 4.6)
(Figure 7D).

Pathway enrichment placed oxidative phosphorylation as the top term in
female-labelled Kenyon cells (adjusted *p* = 6.5 × 10⁻¹⁵), with axon guidance and
cell adhesion in the male-labelled arm (Figure 8). Under the earlier, incorrect
label assignment the same analysis placed oxidative phosphorylation in the male
arm — illustrating that this enrichment tracks specific samples rather than
either biological variable.

![](results/figures/05_volcano_male.png)

![](results/figures/05_volcano_female.png)

![](results/figures/fig4c_de_burden.png)

![](results/figures/fig5c_downsampled_ranking.png)

**Figure 7. Differential expression, reported but not interpreted as a cocaine
effect.** Volcano plots of the nominal cocaine-versus-sucrose contrast in (A) the
male-labelled and (B) the female-labelled samples; dotted lines mark |log₂FC| = 1.
(C) Significant genes per cell-type group and declared sex (mixed-identity
cluster omitted). (D) Significant genes per group after downsampling to 400
cells (mean ± SD, three seeds).

![](results/figures/06_enrichment_cluster_Kenyon-cells_male.png)

![](results/figures/06_enrichment_cluster_Kenyon-cells_female.png)

**Figure 8. Pathway enrichment in Kenyon cells.** Top enriched terms among genes
differentially expressed in (A) male-labelled and (B) female-labelled Kenyon
cells. The red line marks adjusted *p* = 0.05.

### 3.5 Statistical controls

Four controls were applied. Their results are reported here because they
characterise the data, though they inherit the interpretive limitation of §3.3.

**Replicate-level aggregation.** Summing counts per sample gave no genes
significant at *n* = 2 versus 2 (minimum adjusted *p* 0.080 in the male-labelled
arm and 0.052 in the female-labelled arm), as expected. Direction agreement with
the cell-level test was 100% in both arms, with Spearman ρ = 0.941
(*p* = 4.5 × 10⁻¹⁸) and 0.957 (*p* = 1.9 × 10⁻²⁰) over 37 genes each
(Figure 9A, B).

**Interaction model.** Fitting `expression ~ sex + treatment + sex:treatment` to
pseudobulk profiles returned **no gene with a significant interaction term** at
FDR < 0.05 (minimum adjusted *p* = 0.75; Figure 9C). With four residual degrees of
freedom this is expected rather than informative. It does mean that the sexual
dimorphism reported by the original study, which rests on comparing separate
per-sex gene lists rather than on an interaction test, is not demonstrated for
any individual gene in these data (Gelman and Stern 2006).

The main effects are more informative. Twenty-three genes showed a significant
main effect of treatment at FDR < 0.05. The two strongest were `roX2`
(−8.50 log₂, adjusted *p* = 0.0049) and `roX1` (−7.96, adjusted *p* = 0.0052),
with `Yp1` sixth (+4.68, adjusted *p* = 0.018). The main effect of declared sex on
the same two male-specific lncRNAs was close to zero (`roX1` −0.25, `roX2` +0.10
log₂). At replicate level, therefore, the sex markers are explained by the
treatment label and not by the sex label — the same pattern as §3.3, obtained
without cell-level pseudoreplication.

**Permutation control.** Within each declared sex, four samples permit two
independent null contrasts besides the treatment split — the replicate axis and a
diagonal pairing — each using identical cells, test and thresholds (Figure 9D):

| | Treatment | Replicate axis (null) | Diagonal (null) |
|---|---|---|---|
| Male-labelled | 90 | 16 | 44 |
| Female-labelled | 152 | 19 | 33 |

Both treatment contrasts exceed the mean of their two nulls (3.0-fold and
5.9-fold). Given §3.3, this is expected: the treatment split separates the
samples along the axis the sex markers identify, while the null splits mix them.

**Cluster-size matching** is reported in §3.4.

![](results/figures/05b_pseudobulk_concordance_male.png)

![](results/figures/05b_pseudobulk_concordance_female.png)

![](results/figures/fig5a_interaction_volcano.png)

![](results/figures/fig7_permutation_control.png)

**Figure 9. Statistical controls.** Agreement between cell-level Wilcoxon and
pseudobulk log₂ fold changes for the 37 cell-level significant genes in (A) the
male-labelled and (B) the female-labelled arm. (C) Sex × treatment interaction
coefficients from the pseudobulk model; no gene is significant at FDR < 0.05.
(D) Differentially expressed genes for the treatment split compared with two null
splits of the same samples.

### 3.6 Partial reproduction of Baker et al. Figure 3

At the paper's threshold, the gene × cell-type matrices contained 168 genes across
23 cell-type groups in the male-labelled arm and 88 across 23 in the
female-labelled arm, against 133 and 54 reported (Figure 10A). Venn analysis
across six cell types gave Jaccard indices of 0.048–0.116, consistent with the
limited overlap the original describes (Figure 10B).

![](results/figures/fig6ab_gene_cluster_matrix.png)

![](results/figures/fig6c_sex_overlap_venn.png)

**Figure 10. Partial reproduction of Baker et al. Figure 3.** (A) Genes
differentially expressed at the paper's threshold (|log₂FC| > 1.443) in each
cell-type group, in the male-labelled (top) and female-labelled (bottom) arms;
magenta, up in cocaine-labelled samples; cyan, down. (B) Overlap between the two
arms in six cell types matched to the paper's clusters by marker genes.

---

## 4. Discussion

### What reproduces

Three published results reproduce under an independent pipeline. Cell count after
quality control matched to within 0.05%, confirming that the filtering criteria
and input data were correctly identified. Clustering recovered the major neuronal
and glial populations, and the embedding shows no evidence of batch effect — only
one of 30 clusters exceeded 30% composition from a single sample, and that
cluster proved to be a genuine biological population. Cell types corresponding to
the paper's Kenyon cells, surface glia and unannotated antennal/optic cluster
were identified independently by marker gene.

The atlas, in other words, reproduces. What cannot be reproduced from these data
is any conclusion about cocaine.

### The labelling discrepancy

Sex-specific markers place all four sucrose-labelled samples in one group and all
four cocaine-labelled samples in the other. The evidence is not marginal: `roX1`
and `roX2` differ by more than 300-fold between the groups at replicate level and
roughly 700-fold in the cell-level test, and yolk proteins, which males do not
transcribe, read 0.010 or less in every sucrose sample. The replicate-level model
reaches the same conclusion: the treatment label explains `roX1` and `roX2`, and
the declared sex label explains almost none of their variation.

The sample identities themselves are not in doubt. Per-sample cell counts in
Supplemental Table S2 match the eight deposited matrices exactly; the GEO titles
agree; and the authors' Supplemental Code documents the same assignment from
CellRanger output directories onward. The discrepancy therefore did not arise
during analysis or deposition — in either their pipeline or ours.

Neither analysis checked it. The authors' code assigns sex as a hard-coded string
based on which output directory a matrix came from; no marker-based verification
appears anywhere in it. Our own pipeline, as first written, did the same. The
check that revealed the problem — comparing `roX1`/`roX2` and `Yp1`–`Yp3` against
the declared labels — takes under a minute and is not, at present, standard
practice. The signal was nonetheless visible in the original analysis: `roX2`
appears in the published list of genes responding globally to cocaine.

We do not claim to know which explanation is correct. A transposition of the sex
annotation before sequencing would account for it, as would some aspect of the
experimental design not described in the methods. What we can say is that in the
data as deposited, the cocaine-versus-sucrose contrast within each declared sex
separates samples along an axis that sex markers identify as sex. On that basis
we decline to interpret any differential expression result here as a drug effect,
and we suggest the same caution applies to the corresponding results in the
original report until the discrepancy is resolved.

### What this does not show

It does not show that the biological conclusions of Baker et al. are wrong. Acute
cocaine may well produce cell-type-specific and sexually dimorphic transcriptional
change in the fly brain; the behavioural data in that paper, which we did not
reanalyse, are independent of this issue. The claim here is narrower: that the
deposited data cannot be used to demonstrate it without first resolving the
labelling question.

Nor does it show that the original analysis was careless. The pipeline is
documented, the code is public, and the supplementary material is detailed enough
that we were able to verify sample identity three independent ways — which is
more than many published datasets permit, and is the reason this discrepancy was
findable at all.

### Methodological observations

Several observations stand independently of the labelling question.

**Marker verification should be routine.** A dataset with a declared sex variable
can be checked against sex markers in one line of code. The same applies to any
declared variable with a known transcriptional signature. Neither the original
analysis nor ours did this by default.

**Reference gene symbols require checking.** Five fly symbols in this reference
are replaced by vertebrate ortholog names, and one substitution (`trh` =
*trachealess*) silently returns an unrelated gene under case-insensitive
matching. Marker panels taken from the literature will not match a reference
without verification by gene ID.

**Gene counts need a background.** Reporting that a contrast yields *n*
differentially expressed genes is uninterpretable without knowing what an
equivalent incorrect contrast yields. In these data, null contrasts within each
declared sex returned 16 to 44 genes at the same thresholds.

**Cell-level differential expression needs replicate-level support.** Applying a
per-cell test across tens of thousands of cells from two animals per group
inflates significance (Squair et al. 2021; Zimmerman et al. 2021). Pseudobulk
aggregation here returned no significant genes within either declared sex, as
expected at *n* = 2.

### Limitations of this reanalysis

**The discrepancy is unresolved.** We have not determined whether the sex
annotation is wrong or whether something else explains the pattern. That
determination requires information the deposited data does not contain.

**An earlier version of this analysis applied an incorrect label reassignment**,
inferred from the deposition order in the authors' published code. Supplemental
Table S2 showed it to be wrong and it was reverted; the history is preserved in
the repository. It is included here because the error is instructive — a
plausible reconstruction from one source was refuted by another, and only the
cell counts, which are intrinsic to the files, settled it.

**Six clusters (15.7% of cells) are unannotated**, including the cluster ranking
second under size-matched analysis.

**The paper's astrocyte finding could not be assessed.** The genuine
astrocyte-like population identified here contains 143 cells, too few for
differential expression.

**Enrichment used a whole-genome background**, which is optimistic for
brain-expressed terms.

### Conclusion

An independent reanalysis of GSE152495 reproduces the published cell count to
within 0.05% and recovers the cell-type atlas, but finds that sex-specific
markers in the deposited data do not agree with the deposited sex labels. The
consequence is that neither treatment contrast can be separated from a sex
contrast, and no differential expression result from these data — ours or the
original — can presently be attributed to cocaine. The question has been referred
to the original authors, and this report will be updated when they respond.

---

## References

Baker BM, Mokashi SS, Shankar V, Hatfield JS, Hannah RC, Mackay TFC, Anholt RRH. 2021. The *Drosophila* brain on cocaine at single-cell resolution. *Genome Res* 31: 1927–1937. doi:10.1101/gr.268037.120

Becker JB, Koob GF. 2016. Sex differences in animal models: focus on addiction. *Pharmacol Rev* 68: 242–263. doi:10.1124/pr.115.011163

Chang H, Grygoruk A, Brooks E, Ackerson LC, Maidment NT, Bainton RJ, Krantz DE. 2006. Overexpression of the *Drosophila* vesicular monoamine transporter increases motor activity and courtship but decreases the behavioral response to cocaine. *Mol Psychiatry* 11: 99–113. doi:10.1038/sj.mp.4001742

Corey JL, Quick MW, Davidson N, Lester HA, Guastella J. 1994. A cocaine-sensitive *Drosophila* serotonin transporter: cloning, expression, and electrophysiological characterization. *Proc Natl Acad Sci* 91: 1188–1192. doi:10.1073/pnas.91.3.1188

Croset V, Treiber CD, Waddell S. 2018. Cellular diversity in the *Drosophila* midbrain revealed by single-cell transcriptomics. *eLife* 7: e34550. doi:10.7554/eLife.34550

Davie K, Janssens J, Koldere D, De Waegeneer M, Pech U, Kreft Ł, Aibar S, Makhzami S, Christiaens V, Bravo González-Blas C, et al. 2018. A single-cell transcriptome atlas of the aging *Drosophila* brain. *Cell* 174: 982–998.e20. doi:10.1016/j.cell.2018.05.057

Fang Z, Liu X, Peltz G. 2023. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics* 39: btac757. doi:10.1093/bioinformatics/btac757

Filošević A, Al-Samarai S, Andretić Waldowski R. 2018. High throughput measurement of locomotor sensitization to volatilized cocaine in *Drosophila melanogaster*. *Front Mol Neurosci* 11: 25. doi:10.3389/fnmol.2018.00025

Gelman A, Stern H. 2006. The difference between "significant" and "not significant" is not itself statistically significant. *Am Stat* 60: 328–331. doi:10.1198/000313006X152649

Highfill CA, Baker BM, Stevens SD, Anholt RRH, Mackay TFC. 2019. Genetics of cocaine and methamphetamine consumption and preference in *Drosophila melanogaster*. *PLoS Genet* 15: e1007834. doi:10.1371/journal.pgen.1007834

Kuleshov MV, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z, Koplev S, Jenkins SL, Jagodnik KM, Lachmann A, et al. 2016. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Res* 44: W90–W97. doi:10.1093/nar/gkw377

Lukas SE, Sholar M, Lundahl LH, Lamas X, Kouri E, Wines JD, Kragie L, Mendelson JH. 1996. Sex differences in plasma cocaine levels and subjective effects after acute cocaine administration in human volunteers. *Psychopharmacology* 125: 346–354. [CHECK before submission]

McClung C, Hirsh J. 1998. Stereotypic behavioral responses to free-base cocaine and the development of behavioral sensitization in *Drosophila*. *Curr Biol* 8: 109–112. doi:10.1016/s0960-9822(98)70041-7

McInnes L, Healy J, Melville J. 2018. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. arXiv:1802.03426.

Pandey UB, Nichols CD. 2011. Human disease models in *Drosophila melanogaster* and the role of the fly in therapeutic drug discovery. *Pharmacol Rev* 63: 411–436. doi:10.1124/pr.110.003293

Squair JW, Gautier M, Kathe C, Anderson MA, James ND, Hutson TH, Hudelle R, Qaiser T, Matson KJE, Barraud Q, et al. 2021. Confronting false discoveries in single-cell differential expression. *Nat Commun* 12: 5692. doi:10.1038/s41467-021-25960-2

Thurmond J, Goodman JL, Strelets VB, Attrill H, Gramates LS, Marygold SJ, Matthews BB, Millburn G, Antonazzo G, Trovisco V, et al. 2019. FlyBase 2.0: the next generation. *Nucleic Acids Res* 47: D759–D765. doi:10.1093/nar/gky1003

Traag VA, Waltman L, van Eck NJ. 2019. From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep* 9: 5233. doi:10.1038/s41598-019-41695-z

Wang KH, Penmatsa A, Gouaux E. 2015. Neurotransmitter and psychostimulant recognition by the dopamine transporter. *Nature* 521: 322–327. doi:10.1038/nature14431

Wolf FA, Angerer P, Theis FJ. 2018. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19: 15. doi:10.1186/s13059-017-1382-0

Young MD, Behjati S. 2020. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. *GigaScience* 9: giaa151. doi:10.1093/gigascience/giaa151

Zimmerman KD, Espeland MA, Langefeld CD. 2021. A practical solution to pseudoreplication bias in single-cell studies. *Nat Commun* 12: 738. doi:10.1038/s41467-021-21038-1

---

## Supplementary figures

![](results/figures/fig1_qc_scatter.png)

**Supplementary Figure S1.** UMI counts against genes detected per cell before
filtering, coloured by mitochondrial fraction; dashed lines mark the gene-count
thresholds.

![](results/figures/figS_pca_elbow.png)

**Supplementary Figure S2.** Variance explained by the first 50 principal
components; 30 components (16.9% of variance) were retained.

![](results/figures/fig2c_umap_by_sex.png)

![](results/figures/fig2d_umap_by_treatment.png)

**Supplementary Figure S3.** UMAP split by (A) declared sex and (B) treatment.

![](results/figures/figS_composition_sex_treatment.png)

**Supplementary Figure S4.** Composition of each cluster by declared sex (top)
and treatment (bottom); dashed lines mark the dataset-wide proportion.

---

## Appendix A. AI Usage Disclosure

[TO COMPLETE — this appendix must be filled in by the author. The template in
`reports/AI_USAGE_DISCLOSURE.md` requires: each tool and model used and its
purpose; the exact prompts, with date, purpose and what was changed or verified
for each; the completed validation checklist; and the signed statement.]
