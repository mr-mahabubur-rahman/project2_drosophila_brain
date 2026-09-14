# Introduction — draft

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
   model, exclusion of technically labile gene classes, cluster-size matching,
   and comparison against alternative groupings of the same samples; and
4. report clearly which conclusions are supported by these data, which are
   confounded, and which the design cannot resolve.

---

## Notes

- **Check every citation.** I have given them from memory of the paper's
  reference list; verify each against the original before submitting. The ones I
  am least certain of are Squair et al. 2021 (pseudoreplication in single-cell
  DE, *Nature Communications*) and Young & Behjati 2020 (SoupX, *GigaScience*) —
  both are real and relevant, but confirm the details.
- **§1.4 is where this Introduction differs** from a conventional one. It sets
  up verification as a distinct aim from reproduction, which is what makes the
  Results follow naturally rather than appearing as a list of complaints.
- **Consider adding a sentence** on why n = 2 per group was standard practice
  for scRNA-seq in 2020–21. It contextualises the limitation without excusing
  it, and it is fairer to the original authors.
