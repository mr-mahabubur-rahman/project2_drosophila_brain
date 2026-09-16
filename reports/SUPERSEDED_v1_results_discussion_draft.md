# Results and Discussion — draft

*Numbers below come from your own output. Every claim is traceable to a table in
`results/tables/` or a figure in `results/figures/`. Edit freely — you are the
one who has to defend each sentence.*

---

## Results

### 3.1 Three data-handling errors were identified before analysis

Three problems in the supplied data were found and corrected. Each would have
propagated silently; none caused the pipeline to fail.

The supplied `features.tsv.gz` files were space-delimited rather than
tab-delimited, which does not conform to the CellRanger specification that
`scanpy.read_10x_mtx` parses. All three fields collapsed into one column and
gene-symbol lookup failed. Conversion preserved field content and changed only
the delimiter.

The gene *nanchung* (`Dmel_CG5842`, symbol `nan`) is parsed by pandas as a
missing value rather than text, leaving one gene unnamed. It was renamed
`CG5842`, which cannot be misparsed on CSV round-trips. The related gene
*narrow abdomen* (`na`, `Dmel_CG1517`) is unaffected, as lowercase `na` is not
in the default NA list.

Most consequentially, **the folder labels did not match their contents.**
Sex-specific markers separated the eight samples into two groups of four, but
the split followed the treatment field of the folder name rather than the sex
field (Table 1). The male panel (`roX1`, `roX2`, MSL complex) gave a 10.6-fold
non-overlapping separation by treatment label and 1.0-fold by sex label; the
female panel (`Sxl`, `Yp1`–`Yp3`, `tra`) gave 2.1-fold and 1.2-fold
respectively. Yolk proteins read 0.000–0.001 in all four Sucrose-labelled
samples and 0.028–0.422 in all four Cocaine-labelled samples; as yolk proteins
are transcribed only in females, this cannot be a drug effect.

The cause was identified from the deposition order listed in the authors'
published analysis code: pairing that order against the supplied folder names
sorted alphabetically reproduces the observed marker pattern for all eight
samples. Four sex labels and four treatment labels were corrected accordingly.

**Sex assignment is confirmed**: after correction, `roX1` reads 3.898 in males
versus 0.074 in females, and `Yp1`–`Yp3` are ~100-fold higher in females.
**Treatment assignment is inferred** from the deposition order, since no marker
gene reports drug consumption; it is supported independently by the male-biased
response reproducing (§3.4).

### 3.2 Quality control reproduces the published cell count

After filtering (300–2500 genes per cell, ≥5 cells per gene, <10% mitochondrial
reads), **86,177 cells** were retained from 88,991 (96.8%), against the 86,224
reported by Baker et al. — a difference of 47 cells (0.05%). Gene count fell
from 17,481 to 12,189.

Sequencing depth varied roughly two-fold between samples (median 675–1,195 genes
per cell; Figure 1). Critically, depth was **not** confounded with treatment:
within each sex, one deep and one shallow replicate fell in each arm. Had both
deep samples been cocaine, the differential expression results would have been
largely a depth artefact.

### 3.3 Clustering recovers 30 populations; 24 were assignable

Leiden clustering at the published resolution of 0.8 yielded **30 clusters**
against the 36 reported (Figure 2a). A resolution sweep from 0.4 to 1.4 gave a
monotonic increase (19, 24, 30, 33, 34, 37; Figure 2b), reaching 36 at
approximately resolution 1.33. The difference is attributable to method: the
original used Seurat with SCTransform normalization and SNN clustering; this
reanalysis used Scanpy with CPM+log1p and Leiden.

Notably, the published stability argument — that cluster count plateaus at
resolution 0.8 — was not reproduced. Cluster count rose steadily across the
range tested, with no plateau at 0.8.

The embedding is not driven by batch, sex or treatment. Only 1 of 30 clusters
exceeded 30% composition from a single sample (cluster 11, 34.3%, against 12.5%
expected under even mixing). Treatment composition was tight across clusters
(overall 51.8%; maximum deviation cluster 26 at 68.7%). Sex composition varied
more (overall 49.2%; range 22.2% to 71%), identifying a minority of genuinely
sex-biased populations.

Of the 30 clusters, **24 received cell-type labels** and 6 were left
unannotated, comprising 13,526 cells (15.7%). Assignments used two independent
marker panels: the top-three markers published for each of the paper's 36
clusters, and the canonical panel from the project guide (Figure 3b).

Four assignments were verified against both panels:

| Paper cluster | Cell type | This study | Panel score vs. next |
|---|---|---|---|
| C11 | Kenyon cells | cluster 12 | 2.44 vs 2.07 |
| C16 | Unannotated (antennal/optic lobe) | cluster 15 | 2.74 vs 0.96 |
| C22 | Surface glia and fat body | cluster 25 | 2.94 vs 0.96 |

The Kenyon-cell assignment illustrates why two panels were necessary. Ranked on
the published triplet (`jdp`/`Pka-R2`/`Rgk1`) alone, cluster 28 placed a close
second (2.07 vs 2.44) and could plausibly have been assigned as a Kenyon
subpopulation. The independent panel excluded it: `ey` reads 1.084 in cluster 12
and 0.017 in cluster 28 — indistinguishable from background.

The cluster corresponding to the paper's C17 (astrocytes) could not be
identified with confidence. Cluster 24 matched on `Eaat1` and `CG15522` but was
`repo`-negative (0.028) and `alrm`-negative (0.033) with high neuronal markers,
and was left unannotated. The genuine astrocyte-like population (cluster 2:
`repo` 0.830, `alrm` 4.106, `Gs2` 4.829, `Eaat1` 3.642) comprises only 143
cells.

### 3.4 Sexual dimorphism reproduces; the responding genes do not

Males showed **90 differentially expressed genes** and females **14**
(|log₂FC| > 1, BH-adjusted *p* < 0.05), a ratio of 6.4-fold. Baker et al. report
691 versus 322 (2.1-fold). Direction and asymmetry reproduce; absolute counts
differ as expected given Wilcoxon versus MAST and BH versus Bonferroni
correction.

This ratio is, however, the least robust of several that can be computed from
the same data (§3.5), and the female term in it does not survive comparison
against its own background (§3.9). Both qualifications are developed below.

Overlap between sexes was near-zero: 89 male-only genes, 13 female-only, and
1 shared. The original reported limited overlap; this analysis finds less still.

Ranked by raw count, **Kenyon cells responded most strongly among identified
cell types** (96 genes), reproducing the paper's central claim. However, cluster
size correlates with DE count (Spearman ρ = 0.581), reflecting greater
statistical power in larger clusters. Normalised per 1,000 cells, surface glia
and fat body ranks highest among identified types (58.4) while Kenyon cells
falls to seventh (25.1). Both populations are among the paper's four strongest
responders; the original does not appear to adjust for cluster size.

### 3.5 The male bias depends on whether cells are pooled

Summing across the 23 annotated clusters gives **616 DE genes in males and 373
in females, a ratio of 1.65-fold** — close to the 2.15-fold reported by Baker
et al. and four times smaller than the 6.43-fold obtained from the same data
analysed as a single pooled comparison per sex.

The two tests differ in how a spatially restricted signal propagates. Female
responses are concentrated in a few cell types — Kenyon cells (63 genes), the
antennal/optic lobe cluster (58) and dopaminergic neurons (51) — followed by a
steep drop. Male responses are more evenly distributed across clusters (67, 63,
59, 53, 48, 37, 35, 33). A response confined to one cell type is diluted when
all 42,356 female cells are tested together, whereas a response present across
many cell types survives pooling. The pooled comparison therefore favours
brain-wide responses and understates cell-type-restricted ones.

This was tested directly. Excluding all mitochondrial, respiratory-chain and
ribosomal genes (280 genes) changed the pooled ratio from 6.43-fold to
6.83-fold, with male counts falling only from 90 to 82 (Figure 5b). The
mitochondrial genes are therefore the most statistically significant male hits
but are **not** the cause of the count asymmetry, and no technical explanation
for the pooled–per-cluster discrepancy is supported by the data.

The bias also reverses in specific cell types. Females exceeded males in **6 of
23 clusters**, including several of the paper's strongest responders: Kenyon
cells (63 versus 33), the antennal/optic lobe cluster corresponding to their C16
(58 versus 17), dopaminergic neurons (51 versus 35) and surface glia and fat
body (34 versus 14). The aggregate male bias therefore does not describe every
population.

### 3.6 The leading genes in each sex coincide with technical covariates

**In males**, the most significant genes were dominated by mitochondrial and
mitochondria-associated transcripts: `mt:lrRNA` (+3.05 log₂FC), `Cyt-c-p`
(cytochrome c), `Mic10b` (MICOS inner-membrane complex), `sun` (*stunted*, an
ATP synthase subunit) and `VhaM9.7-a` (V-ATPase subunit), alongside
`RNaseMRP:RNA`, `snRNA:7SKa` and `RpL41`. Four nuclear-encoded mitochondrial
genes and the mitochondrial rRNA therefore rise together.

Mitochondrial content differs systematically between treatment arms in males:
median 0.373% in cocaine versus 0.058% in sucrose, a 6.4-fold difference
(Figure 1, panel 3). Females show no such difference (0.260% versus 0.267%) and
no `mt:lrRNA` effect (−0.03).

`RpL41` is an exception worth noting: it is named by Baker et al. among their
globally cocaine-responsive genes and replicates here.

**The scope of this confound is limited.** Excluding all 280 mitochondrial,
respiratory-chain and ribosomal genes reduced the male count from 90 to 82 and
left the male:female ratio unchanged (6.43-fold to 6.83-fold). The mitochondrial
signal therefore affects which genes rank highest, and the interpretation of any
individual mitochondrial gene or of the oxidative-phosphorylation enrichment
below, but it does not account for the magnitude of the male response.

Pathway enrichment reflects the same pattern. Oxidative phosphorylation and
mitochondrial electron transport were the top enriched terms in three separate
male clusters — dopaminergic neurons (adjusted *p* = 5.4 × 10⁻⁸), Kenyon cells,
and neuropeptide/cholinergic neurons. No respiration terms reach significance in
any female cluster; female enrichments were axon guidance, synaptic target
recognition and phototransduction. One respiratory subunit (`UQCR-6.4`) does
appear among the top female genes, so the contrast is at the pathway level
rather than absolute.

**In females**, the top genes were the yolk proteins `Yp1`–`Yp3`. These are
ambient RNA, not cellular expression. The detection rate was 27.2% of cells in
`Female_Cocaine_R2` against 3.2–4.9% in the sucrose samples, with
`Female_Cocaine_R1` intermediate at 7.0%. The signal is distributed across all
cell types, including photoreceptors (`Yp3` mean 0.663) and Kenyon cells
(0.226) — populations that do not transcribe yolk proteins. This is the
signature of free transcripts from lysed fat body partitioning into droplets,
with one dissection disproportionately affected.

`ninaE` (rhodopsin) was also differentially expressed (+1.29 in males), in the
opposite direction to the published result. As a photoreceptor gene, its
abundance in a brain dissociation depends on retinal carry-over rather than
treatment.

### 3.7 No individual gene shows a significant sex x treatment interaction

Dimorphism is established above, and in the original study, by comparing gene
lists derived from separate per-sex analyses. That comparison does not test
whether any gene responds *differently* by sex: a gene significant in one
analysis and not the other may not differ significantly between them (Gelman &
Stern 2006).

Fitting `expression ~ sex + treatment + sex:treatment` to pseudobulk profiles
(9,397 genes, 8 samples) returned **no genes with a significant interaction
term** at FDR < 0.05, and none with a significant main effect of treatment. Of
the 90 genes called significant in males by separate per-sex testing, none
showed a significant interaction.

With four residual degrees of freedom this is an expected outcome rather than
evidence of absence. Genes with the largest interaction coefficients —
`Mic10b` (2.57), `sun` (1.19), `Lztr1` (1.02) — are plausible candidates but
are not distinguishable from noise at this replicate number. The finding is
that **the design supports dimorphism in aggregate but cannot establish it for
any individual gene**, and the same limitation applies to the original analysis.

### 3.8 The response ranking is unstable at matched statistical power

Downsampling every cluster to 400 cells and repeating over three random seeds
gives a ranking free of cluster-size effects. Surface glia and fat body ranked
third (11.7 ± 2.6 genes), consistent with its third-place position under
per-1000-cell normalisation. Kenyon cells fell to fifth (7.3 ± 2.6) from first
in the raw count.

Standard deviations were large relative to the means — the highest-ranked
cluster varied between 10 and 28 genes across seeds — so ordering beyond the
first few positions is not stable. Surface glia is the only cell type ranking
highly under all three approaches.

### 3.9 The female response is below its own background

A differentially expressed gene count is uninterpretable without a baseline.
Within each sex there are four samples and six possible 2-versus-2 splits. One
is the true treatment split, one its complement, and four cut across the
treatment assignment — these last are nulls, since no treatment difference
exists between the groups they define.

| | Real split | Null splits | Ratio |
|---|---|---|---|
| Male | 90 | 7, 7, 14, 14 (median 10) | **8.6x** |
| Female | 14 | 35, 35, 39, 39 (median 37) | **0.38x** |

The male response exceeds every null by a wide margin: it is specific to
treatment. **The female response falls below every null.** Any alternative
grouping of the same four female samples yields more differentially expressed
genes than the correct one.

Gene identity confirms this. Of the 90 male genes, 7 (8%) also appear in a null
split, and those seven are structural and stress transcripts —
`28SrRNA-Psi:CR40596`, `28SrRNA-Psi:CR45848`, `snRNA:7SKa`, `RNaseMRP:RNA`,
`Hsp26`, `Hsp27`, `Rh2`. Of the 14 female genes, **9 (64%) appear in a null
split**, including `Yp1`, `Yp2` and `Yp3` — the yolk proteins independently
identified as ambient RNA in §3.6, together with `Arr1`, `Cp110`, `Imp`,
`Nplp2` and `RpL39`.

The female result is therefore not a smaller response than the male one. At this
threshold and replicate number it is not distinguishable from variation between
female samples, and the majority of its genes are identifiable as such directly.

### 3.10 Cell-level results are directionally corroborated at replicate level

Pseudobulk analysis (summing counts per sample, n = 2 versus 2) found no
significant genes, as expected at this replicate number. Direction agreement
with the cell-level test was **100% in both sexes**, with Spearman ρ = 0.776
(males, *p* = 3.0 × 10⁻³) and 0.617 (females). The comparison was restricted to
highly variable genes.

---

## Discussion

### The reanalysis reproduces the design-level findings

Three published results reproduce. Cell count after QC matched to within 0.05%,
confirming that the filtering criteria and input data were correctly identified.
The male bias in transcriptional response reproduced at 6.4-fold against a
reported 2.1-fold, with near-zero overlap between sexes. Kenyon cells emerged as
the most responsive identified cell type, matched to the published cluster by
marker gene rather than cluster number.

These are the claims that survive independent reanalysis, and they are the
claims the paper's title and abstract rest on: that the *Drosophila* brain
responds to acute cocaine in a cell-type-specific and profoundly sexually
dimorphic manner.

### The female response is not established by these data

The permutation control is the most consequential result of this reanalysis.
Fourteen differentially expressed genes in females sounds like a modest but real
effect until it is compared with the 35 to 39 genes obtained from any incorrect
grouping of the same four samples. A signal smaller than its own background
cannot be attributed to the variable of interest.

This does not show that female flies do not respond to cocaine. A real effect
smaller than the between-sample variation would look identical, and with two
replicates per arm that variation is substantial. What it shows is that **this
design cannot resolve a female response at this threshold**, and that the
reported female gene count therefore carries less information than its
appearance suggests.

It also explains two earlier results. The 6.43-fold pooled ratio compares a real
male signal against a female floor rather than two responses of different
magnitude. And no gene showed a significant sex x treatment interaction (§3.7)
because the female arm contributes no reliable signal for the male arm to
interact with.

### The gene-level findings are confounded in both sexes

The genes driving those results are a different matter. In males, the most
significant transcripts are mitochondrial, and mitochondrial content differs
6.4-fold between treatment arms. In females, the most significant transcripts
are yolk proteins present as ambient RNA, concentrated in one dissection.

For the female yolk-protein result these are not competing explanations to be
weighed against a biological one — they are the more parsimonious account.
Cocaine does not induce yolk protein transcription in photoreceptors.

The male mitochondrial case is narrower than it first appears. A 6.4-fold
difference in mitochondrial fraction between treatment groups, present in one
sex and absent in the other, is more readily explained by dissociation stress
than by drug action. But excluding mitochondrial and respiratory genes
altogether left the male response count and the male:female ratio essentially
unchanged, so the confound affects the identity and ranking of the top male
genes — and the oxidative-phosphorylation pathway result — without accounting
for the response as a whole.

The female data functions as an internal control throughout. Where males show a
mitochondrial difference between arms, females do not; where males show
respiration pathways enriched in three clusters, females show none. The
asymmetry is consistent with a technical origin in the male dissections.

Neither check appears in the original report. Per-cluster and per-arm
mitochondrial fraction is not tabulated, and ambient RNA contamination is not
assessed. Both are standard quality-control steps in current practice, and both
are computable from the deposited data.

### What this implies for the published conclusions

The conclusion that responses are sexually dimorphic is not threatened; if
anything it is strengthened, since it reproduces under a different pipeline.
The conclusion that specific cell types respond differentially also holds, with
the caveat below about cluster size.

What requires qualification is the gene-level and pathway-level interpretation.
The paper reports `mt:lrRNA` among its globally cocaine-responsive genes and
identifies metabolic pathway remodelling. This reanalysis suggests that at least
part of that signal reflects mitochondrial content differing between the groups
being compared, rather than regulated expression.

### The male bias is real, but its magnitude is analysis-dependent

The pooled and per-cluster ratios differ four-fold (6.43 versus 1.65) on the
same corrected data. A technical explanation was tested and rejected: excluding
mitochondrial, respiratory and ribosomal genes did not reduce the pooled ratio.

The more likely explanation is biological. Female responses concentrate in a
few cell types while male responses distribute across many. Pooling dilutes a
cell-type-restricted response and preserves a distributed one, so the pooled
comparison systematically favours the male pattern. On this reading the sexes
differ not only in how *much* of the brain responds but in how *focally* —
which is a stronger and more specific claim than "males respond more", and one
the original analysis, being pooled, could not have detected.

Reporting both ratios is therefore necessary. Either alone is misleading about
what is being measured.

### Cluster size confounds the response ranking

The correlation between cluster size and DE gene count (ρ = 0.581) means a raw
ranking partly measures statistical power. Surface glia moves from twelfth to
third when normalised; Kenyon cells falls from first to seventh. Both are
defensible rankings of different quantities, and reporting only one risks
presenting power as biology. Cluster-size normalisation, or a fixed-power
approach such as downsampling to equal cell numbers, would strengthen
conclusions of this kind.

### Limitations of this reanalysis

**Pseudoreplication.** The cell-level Wilcoxon test treats ~86,000 cells as
independent replicates when there are two biological replicates per group. The
pseudobulk cross-check showed 100% direction agreement, supporting the
directional findings, but it cannot establish significance at n = 2 and was
restricted to highly variable genes. This limitation applies equally to the
original analysis, which used MAST on individual cells.

**Cluster count.** Thirty clusters were recovered rather than 36, and six
clusters (15.7% of cells) could not be assigned with confidence. Finer
annotation would require reference-based label transfer against an existing fly
brain atlas rather than manual marker inspection.

**Astrocytes.** The published C17 astrocyte response could not be assessed. The
genuine astrocyte-like population identified here contains 143 cells — too few
for differential expression. This is a power limitation, not a contradiction.

**Treatment assignment.** Sex labels were confirmed against marker genes;
treatment labels rest on the deposition order in the authors' repository. No
computational check can verify which animals consumed cocaine.

**Enrichment background.** Enrichr tests against all annotated fly genes,
whereas the differential expression test could only detect brain-expressed
genes surviving QC. Adjusted *p*-values are therefore optimistic for
brain-expressed terms.

### Concluding remark

Independent reanalysis of a published dataset serves two purposes: to test
whether results reproduce, and to test whether they mean what was claimed. Here
the design-level findings reproduce robustly, while the gene-level findings
appear substantially confounded by technical variation that was not reported.
Both outcomes are informative, and neither is visible without recomputing the
analysis from the deposited data.

---

## Notes for you before submitting

- **Check `mt:srRNA` and `RNaseMRP:RNA`** against your volcano and DE tables;
  I read those from the figure labels.
- **The astrocyte paragraph is the weakest point.** A reviewer may ask whether
  cluster 24 is astrocytes after all and your `repo` threshold is too strict.
  Have the numbers to hand: cluster 2 `repo` 0.830 vs cluster 24 at 0.028.
- **Soften if you prefer.** "Suggests that at least part of that signal
  reflects" is the strongest claim I have made; "is consistent with" is the
  weaker alternative if you would rather not defend the stronger one.
- **Cite** Baker et al. 2021 (*Genome Res* 31:1927–1937), Wolf et al. 2018
  (Scanpy), Traag et al. 2019 (Leiden), McInnes et al. 2018 (UMAP), and FlyBase
  for every marker gene asserted.
