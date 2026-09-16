# Rewritten report — Abstract, Methods §2.1, Results, Discussion

*Replaces the corresponding sections of the previous draft. Every number is from
your current output under the deposited sample labels.*

**Framing note:** this is written as an **unresolved discrepancy referred to the
authors**, not as a definitive claim that the published design is confounded.
That is the defensible position until they reply. If they confirm a label
transposition, §3.3 and the Discussion tighten considerably; if they confirm the
labels as correct, the finding becomes stronger, not weaker.

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
`roX2` read 3.8–4.0 and 2.1–2.4 in all four sucrose samples and 0.06–0.11 and
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

---

## Methods §2.1 — Data and sample identity

Count matrices for GSE152495 were obtained pre-processed as eight CellRanger v3
output directories. Three problems were corrected before analysis; each would
have propagated silently.

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
the repository history.

The sex markers nonetheless disagree with these labels (§3.3), a discrepancy
that could not be resolved from the deposited data.

---

## 3. Results

### 3.1 Quality control reproduces the published cell count

After filtering (300–2,500 genes per cell, <10% mitochondrial reads, ≥5 cells per
gene), **86,177 of 88,991 cells were retained (96.8%)**, against the 86,224
reported — a difference of 47 cells, or 0.05%. Gene count fell from 17,481 to
12,189. The agreement indicates the filtering criteria and input data were
correctly identified.

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

Group means suggest a female-specific treatment difference (0.571 versus 0.103)
and a male difference in the opposite direction (0.348 versus 0.528), but both
are driven by single replicates — `Female_Cocaine_R2` and `Male_Sucrose_R2`
respectively. This is dissection-to-dissection variation, not a group-level
property, and any mitochondrial gene appearing in a treatment comparison should
be interpreted with that in mind.

### 3.2 Clustering recovers 30 populations; 24 were assignable

Leiden clustering at resolution 0.8 yielded 30 clusters against the 36 reported.
A resolution sweep gave a monotonic increase (19, 24, 30, 33, 34, 37 at
resolutions 0.4 to 1.4), reaching 36 at approximately 1.33. **The published
stability argument — that cluster count plateaus at resolution 0.8 — was not
reproduced.** This may be a property of Leiden rather than a disagreement about
the data.

The embedding is not driven by batch. Only one of 30 clusters exceeded 30%
composition from a single sample (34.3%, against 12.5% under even mixing), and
that cluster remained 74.5% female among depth-matched cells against a 56.4%
background, indicating a genuine population rather than an artefact.

Twenty-four clusters received cell-type labels using two independent marker
panels — the paper's published top-three markers per cluster, and a canonical
*Drosophila* neuronal and glial panel. Six clusters (13,526 cells, 15.7%) were
left unannotated where the panels disagreed. Populations corresponding to the
paper's Kenyon cells (C11), surface glia (C22) and unannotated antennal/optic
lobe cluster (C16) were identified by marker gene rather than cluster number. The
paper's C17 (astrocytes) could not be confidently identified: the cluster matching
on `Eaat1` was `repo`- and `alrm`-negative, while the genuine astrocyte-like
population contains only 143 cells.

### 3.3 Sex markers do not agree with the deposited sample labels

This is the principal finding and it precludes interpretation of the
differential expression results below.

Three sex-specific markers separate the eight samples cleanly, but the split
follows the **treatment** field of the sample name rather than the **sex** field:

| Marker | All four Sucrose samples | All four Cocaine samples |
|---|---|---|
| `roX1` (male-specific lncRNA) | 3.838 – 3.959 | 0.059 – 0.112 |
| `roX2` (male-specific lncRNA) | 2.149 – 2.385 | 0.008 – 0.015 |
| `Yp1`–`Yp3` (female yolk proteins) | 0.000 – 0.010 | 0.028 – 0.428 |
| `Sxl` (female) | 0.695 – 0.984 | 1.227 – 1.762 |

`roX1` and `roX2` are lncRNAs of the dosage-compensation complex, expressed
almost exclusively in males. Yolk proteins are transcribed in female fat body and
are absent from males. The highest yolk-protein value in the dataset (`Yp3` at
0.428) occurs in a sample labelled male.

Genes whose expression is post-transcriptionally regulated and therefore
present in both sexes at the RNA level — `msl-2`, `mle`, `tra` — show no such
separation, as expected. Only the markers that are genuinely sex-specific at the
transcript level split the samples.

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
females.

Sample identity was verified against Supplemental Table S2, the GEO titles and
the authors' Supplemental Code (§2.1); all three agree, and no computational step
in either their pipeline or ours could have introduced the discrepancy. Neither
their code nor ours, as originally written, checked sex markers against the
sample labels.

Two explanations are consistent with the evidence and cannot be distinguished
from the deposited data: the sex annotation may have been transposed at some
point before sequencing, or sex and treatment may covary for a reason not
described in the methods. The question has been referred to the original authors.

### 3.4 Differential expression results, reported but not interpreted

The results below are recorded for completeness and as evidence bearing on §3.3.
**None should be read as a cocaine effect**, since each contrast is confounded
with sex.

Pooled within each declared sex, the male-labelled contrast returned **90
differentially expressed genes** (75 up, 15 down) and the female-labelled
contrast **152** (126 up, 26 down) at |log₂FC| > 1 and BH-adjusted *p* < 0.05 —
a ratio of 0.59, favouring the female-labelled arm. At the paper's stricter
threshold (|log₂FC| > 1.443) the counts were 20 and 61. Baker et al. report 691
and 322, a ratio of 2.15 favouring males. Summing across the 23 annotated
clusters gives 1,282 and 1,256, a ratio of 1.02.

Overlap between the two arms was 53 genes unique to one, 115 to the other, and 37
shared, of which one changed direction. Per cluster, Jaccard indices ranged from
0.048 to 0.116.

Cluster size correlated with differentially expressed gene count (Spearman
ρ = 0.634), so the raw ranking partly measures statistical power. Downsampling
every cluster to 400 cells over three seeds gave central brain B cholinergic
neurons highest (69.0 ± 8.0), followed by the unannotated Eaat1-high cluster
(43.7 ± 11.3), the antennal/optic lobe cluster (39.7 ± 15.7) and surface glia
(38.3 ± 7.4); Kenyon cells ranked sixth (26.3 ± 4.6).

Pathway enrichment placed oxidative phosphorylation as the top term in
female-labelled Kenyon cells (adjusted *p* = 6.5 × 10⁻¹⁵), with axon guidance and
cell adhesion in the male-labelled arm. Under the earlier, incorrect label
assignment the same analysis placed oxidative phosphorylation in the male arm —
illustrating that this enrichment tracks specific samples rather than either
biological variable.

### 3.5 Statistical controls

Four controls were applied. Their results are reported here because they
characterise the data, though they inherit the interpretive limitation of §3.3.

**Replicate-level aggregation.** Summing counts per sample gave no genes
significant at *n* = 2 versus 2, as expected. Direction agreement with the
cell-level test was 100% in both arms, with Spearman ρ = 0.941 (*p* = 4.5 × 10⁻¹⁸)
and 0.957 (*p* = 1.9 × 10⁻²⁰) over 37 genes each.

**Interaction model.** Fitting `expression ~ sex + treatment + sex:treatment` to
pseudobulk profiles returned **no genes with a significant interaction term** at
FDR < 0.05, and none with a significant main effect of treatment. With four
residual degrees of freedom this is expected rather than informative. It does
mean that the sexual dimorphism reported by the original study, which rests on
comparing separate per-sex gene lists rather than on an interaction test, is not
demonstrated for any individual gene in these data (Gelman and Stern 2006).

**Permutation control.** Within each declared sex, four samples permit two
independent null contrasts besides the treatment split — the replicate axis and a
diagonal pairing — each using identical cells, test and thresholds:

| | Treatment | Replicate axis (null) | Diagonal (null) |
|---|---|---|---|
| Male-labelled | 90 | 16 | 44 |
| Female-labelled | 152 | 19 | 33 |

Both treatment contrasts exceed their nulls (3.0-fold and 5.9-fold). Given §3.3,
this is expected: the treatment split separates the samples along the axis the
sex markers identify, while the null splits mix them.

**Cluster-size matching** is reported in §3.4.

### 3.6 Partial reproduction of Baker et al. Figure 3

At the paper's threshold, the gene × cluster matrices contained 46 genes across
23 cell types in the male-labelled arm and 45 across 19 in the female-labelled
arm, against 133 and 54 reported. Venn analysis across six cell types gave
Jaccard indices of 0.048–0.116, consistent with the limited overlap the original
describes.

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
and `roX2` differ by roughly 700-fold between the groups, and yolk proteins,
which males do not transcribe, read 0.000–0.001 in every sucrose sample.

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
practice.

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
aggregation here returned no significant genes in either arm, as it must at
*n* = 2.

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

## Notes for you

- **Figure numbering** will need redoing; several figures from the old draft no
  longer correspond to any claim. The permutation figure, QC figure and marker
  dotplot all stand. The ratio comparison and gene-exclusion figures relate to
  arguments that no longer appear.
- **A new figure would help §3.3** — a simple per-sample bar chart of `roX1` and
  `Yp1` with the sucrose/cocaine split marked. That is the whole finding in one
  panel. Say the word and I will write it.
- **§3.4 and §3.5 could be moved to an appendix.** They are reported for
  completeness but no longer carry the argument. Keeping them in the main text
  risks a reader treating the numbers as results.
- **Update the graphical abstract and flowchart**, both of which still encode the
  old conclusion.
