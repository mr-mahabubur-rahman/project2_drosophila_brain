# Data provenance — GSE152495

Record of where the data came from, what was wrong with it, what was corrected,
and what could not be resolved. Written so that anyone re-running this analysis
can verify each step independently.

**Status:** the sample labelling discrepancy described in §4 is **unresolved** and
has been referred to the original authors. Until it is resolved, no differential
expression result from this dataset should be interpreted as a cocaine effect.

---

## 1. Source

| | |
|---|---|
| Accession | [GSE152495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152495) |
| Publication | Baker BM, Mokashi SS, Shankar V, Hatfield JS, Hannah RC, Mackay TFC, Anholt RRH (2021) *Genome Research* 31:1927–1937. [doi:10.1101/gr.268037.120](https://doi.org/10.1101/gr.268037.120) |
| Design | 8 samples — sex (F/M) × treatment (cocaine/sucrose) × 2 biological replicates |
| Platform | 10x Genomics Chromium; CellRanger v3.1; *D. melanogaster* Release 6 |
| Obtained | Pre-processed CellRanger output directories (barcodes / features / matrix) |

Count matrices are not redistributed in this repository. `tools/download_data.sh`
and `tools/organize_data.sh` rebuild `data/` from GEO.

---

## 2. Corrections applied before analysis

Three problems were found in the supplied files. Each would have propagated
silently; none caused the pipeline to fail.

### 2.1 Space-delimited `features.tsv.gz`

All eight `features.tsv.gz` files were space- rather than tab-delimited. The
CellRanger specification that `scanpy.read_10x_mtx` parses expects tabs, so all
three fields collapsed into a single column and gene-symbol assignment failed.

**Fix:** `tools/fix_features_separator.py --apply` converts the delimiter with
field content unchanged (17,481 lines per file). Originals are preserved in
`data/<sample>/original/`.

### 2.2 Gene `nan` parsed as a missing value

*nanchung* (Dmel_CG5842, FlyBase symbol `nan`) is one of pandas' default
missing-value strings. On read it became `NaN`, leaving one gene unnamed and
preventing HDF5 serialisation of the AnnData object.

**Fix:** renamed to `CG5842` in `config.GENE_NAME_NA_FIXES`. The identifier
survives CSV round-trips, unlike `nan`.

*narrow abdomen* (`na`, Dmel_CG1517) is unaffected — lowercase `na` is not in the
default list.

### 2.3 Vertebrate ortholog names in the reference

The CellRanger reference used for this dataset substitutes vertebrate ortholog
names for five fly symbols. Resolved by gene ID against FlyBase (Thurmond et al.
2019):

| Reference symbol | Actual fly gene | Gene ID |
|---|---|---|
| `VGlut1` | VGlut | CG9887 |
| `Adcy1` | *rut* (rutabaga) | CG9533 |
| `Pde4` | *dnc* (dunce) | CG32498 |
| `Trhn` | Trh (tryptophan hydroxylase) | CG9122 |
| `Tret1` | Tret1-1 | CG30035 |

**One of these is a trap.** `trh` in lowercase is **Dmel_CG42865,
*trachealess*** — a tracheal transcription factor unrelated to serotonin
synthesis. A case-insensitive marker match returns *trachealess* and would label
a cluster serotonergic on the strength of a tracheal gene.

**Fix:** `config.REFERENCE_SYMBOL_ALIASES`, with case-sensitive matching
throughout. `tools/resolve_markers.py` checks every marker in the panel against
the actual annotation before use: 31 matched directly, 5 via alias, 0 dangerous
matches, 0 missing.

---

## 3. Sample identity — verified three ways

Sample assignment was checked against three independent sources. **All three
agree.**

### 3.1 Supplemental Table S2 cell counts

Table S2 of the paper lists per-sample cell counts. These are intrinsic to each
matrix file, so they tie every folder to the sample name the authors used.

| Table S2 label | Cells | Folder as supplied | Cells | Match |
|---|---|---|---|---|
| ♀ Sucrose Rep 1 | 9,072 | Female_Sucrose_1 | 9,072 | ✓ |
| ♀ Sucrose Rep 2 | 11,693 | Female_Sucrose_2 | 11,693 | ✓ |
| ♂ Sucrose Rep 1 | 13,193 | Male_Sucrose_1 | 13,193 | ✓ |
| ♂ Sucrose Rep 2 | 11,033 | Male_Sucrose_2 | 11,033 | ✓ |
| ♀ Cocaine Rep 1 | 13,072 | Female_Cocaine_1 | 13,072 | ✓ |
| ♀ Cocaine Rep 2 | 9,367 | Female_Cocaine_2 | 9,367 | ✓ |
| ♂ Cocaine Rep 1 | 10,437 | Male_Cocaine_1 | 10,437 | ✓ |
| ♂ Cocaine Rep 2 | 11,124 | Male_Cocaine_2 | 11,124 | ✓ |
| **Total** | **88,991** | | **88,991** | ✓ |

Eight of eight, and the total matches the merged dataset exactly.

*(Table S2 contains a typographical error: the fourth row is labelled
"♂ Sucrose Rep 1" where "Rep 2" is intended. The cell count is unambiguous.)*

### 3.2 GEO sample titles

The GEO sample titles correspond to the same assignment.

### 3.3 The authors' Supplemental Code

The published Supplemental Code documents the mapping from CellRanger output
directories to Seurat objects:

```
S1   → Female_sucrose_R1        S5   → Female_cocaine_R1
S2   → Female_sucrose_R2        S6   → Female_cocaine_R2
S3   → Male_sucrose_R1          S7   → Male_cocaine_R1
S4   → Male_sucrose_R2          S8v2 → Male_cocaine_R2
```

Sex and treatment are applied there as hard-coded strings
(`$gender_stim <- "male_cocaine"`) based on which directory a matrix came from.
No marker-based verification appears anywhere in that code.

**Conclusion:** the labels are used as deposited. The labelling chain is
consistent from CellRanger output directory through to the matrices in this
repository.

---

## 4. Unresolved: sex markers do not agree with the sample labels

### 4.1 The observation

Three sex-specific markers separate the eight samples cleanly, but the split
follows the **treatment** field of the sample name rather than the **sex** field.
Mean log-normalised expression per sample, post-QC:

| Marker | Specificity | 4 Sucrose samples | 4 Cocaine samples |
|---|---|---|---|
| `roX1` | male (dosage compensation lncRNA) | 3.846 – 3.972 | 0.059 – 0.116 |
| `roX2` | male (dosage compensation lncRNA) | 2.160 – 2.396 | 0.008 – 0.015 |
| `Yp1` | female (yolk protein) | 0.000 – 0.001 | 0.065 – 0.255 |
| `Yp2` | female (yolk protein) | 0.003 – 0.003 | 0.028 – 0.153 |
| `Yp3` | female (yolk protein) | 0.006 – 0.010 | 0.069 – 0.428 |
| `Sxl` | female | 0.695 – 0.984 | 1.227 – 1.762 |

The ranges do not overlap on any marker. No sample bridges the gap.

Per sample:

| Sample | roX1 | roX2 | Yp1 | Yp2 | Yp3 |
|---|---|---|---|---|---|
| Female_Sucrose_R1 | 3.968 | 2.396 | 0.001 | 0.003 | 0.006 |
| Female_Sucrose_R2 | 3.846 | 2.160 | 0.000 | 0.003 | 0.008 |
| Female_Cocaine_R1 | 0.062 | 0.008 | 0.066 | 0.028 | 0.069 |
| Female_Cocaine_R2 | 0.116 | 0.015 | 0.065 | 0.056 | 0.091 |
| Male_Sucrose_R1 | 3.850 | 2.310 | 0.001 | 0.003 | 0.010 |
| Male_Sucrose_R2 | 3.972 | 2.368 | 0.001 | 0.003 | 0.008 |
| Male_Cocaine_R1 | 0.059 | 0.011 | 0.104 | 0.062 | 0.161 |
| Male_Cocaine_R2 | 0.068 | 0.013 | 0.255 | 0.153 | 0.428 |

Yolk proteins are transcribed in female fat body and are absent from males. The
highest yolk-protein value in the dataset (`Yp3` = 0.428) is in a sample labelled
male.

### 4.2 Controls that behave as expected

Genes regulated post-transcriptionally, and therefore present in both sexes at
the RNA level, show no such separation: `msl-2` (0.039–0.080), `mle`
(0.252–0.376), `tra` (0.045–0.102). Only markers that are genuinely sex-specific
at the transcript level split the samples, which is the expected behaviour and
argues against a generic technical artefact.

### 4.3 The consequence

Comparing cocaine with sucrose *within each declared sex* returns the sex markers
among the most significant genes in **both** comparisons:

| Gene | Male-labelled contrast | Female-labelled contrast |
|---|---|---|
| `roX2` | −9.620 | −9.629 |
| `roX1` | −9.539 | −9.125 |
| `Yp1` | +7.684 | +6.453 |
| `Yp3` | +5.234 | +3.479 |
| `Yp2` | +5.217 | +3.829 |
| `Sxl` | +1.065 | +1.588 |

log₂ fold changes, all with adjusted *p* at or below machine precision.

Each nominal treatment contrast is separating males from females.

### 4.4 What this does and does not establish

**Established:** the sex markers in the deposited data do not agree with the
deposited sex labels, and the discrepancy did not arise during analysis or
deposition — the labelling chain is consistent from CellRanger output through to
GEO, and was verified three independent ways.

**Not established:** which explanation is correct. Two are consistent with the
evidence and cannot be distinguished from the deposited data:

1. the sex annotation was transposed at some point before sequencing; or
2. sex and treatment covary for a reason not described in the methods.

**Referred to:** the corresponding authors (R.R.H. Anholt, T.F.C. Mackay),
Clemson University.

### 4.5 Reproducing this check

```bash
python tools/make_sexmarker_figure.py
```

Prints the per-sample table above and writes
`results/figures/fig_sexmarker_discrepancy.png`.

---

## 5. Record of an incorrect correction, and its reversion

Documented because the error is instructive.

**What happened.** On finding the marker discrepancy (§4), an earlier version of
this analysis concluded that the folder labels were wrong and reconstructed the
"true" assignment from the deposition order in the authors' published R code —
pairing S1–S8 against the folder names sorted alphabetically. That reconstruction
produced a mapping under which the markers agreed with the sex labels, which
appeared to confirm it.

**Why it was wrong.** Supplemental Table S2 lists per-sample cell counts. Under
the reconstructed mapping, `Female_Cocaine_1` (13,072 cells) would have been
♀ Sucrose Rep 1 — but Table S2 gives that sample 9,072 cells. All eight cell
counts contradicted the reconstruction and matched the folders as supplied.

**What changed as a result.** The reverted analysis gives different numbers
throughout:

| | Incorrect mapping | Deposited labels |
|---|---|---|
| Pooled male:female DE ratio | 6.43 | 0.59 |
| Female arm vs its permutation background | 0.38× (below) | 5.85× (above) |
| Pseudobulk concordance (ρ) | 0.78 / 0.62 | 0.94 / 0.96 |
| Mitochondrial signal | male-specific | female-specific, traced to single samples |

**The lesson.** A plausible reconstruction from one source was refuted by
another. Cell count settled it because it is intrinsic to the file, whereas
deposition order is an inference about how files were ordered. When two sources
disagree, prefer the one that cannot have been reordered.

**Reversion:** `tools/revert_remap.py --apply`. The previous configuration is
preserved at `scripts/config.py.pre-revert`, and outputs from the earlier
analysis are marked `SUPERSEDED_` in `reports/`.

---

## 6. What reproduces

Independent of the labelling question:

| | Published | This analysis |
|---|---|---|
| Cells after QC | 86,224 | 86,177 (0.05% difference) |
| Genes retained | — | 12,189 of 17,481 |
| Clusters at resolution 0.8 | 36 | 30 (36 at ≈1.33) |
| Cell types assigned | 36 | 24 of 30; 6 unannotated (15.7% of cells) |
| Batch effect | none reported | 1 of 30 clusters >30% from one sample |

Clustering and annotation do not depend on sample labels and are unaffected by
§4.

---

## 7. Recommendations for anyone using this dataset

1. **Check sex markers against declared labels before analysis.** One line of
   code; neither the original analysis nor this one did it by default.
2. **Verify marker gene symbols by gene ID**, not by name. This reference
   contains at least one substitution (`trh` = *trachealess*) that silently
   returns an unrelated gene.
3. **Verify sample identity by cell count** against Supplemental Table S2 before
   accepting any relabelling.
4. **Do not interpret differential expression from these data as a cocaine
   effect** until the labelling question is resolved.
