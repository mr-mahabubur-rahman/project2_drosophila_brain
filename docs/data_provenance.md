# Data provenance

Record of the data's origin and of three corrections applied before analysis.
All three were silent failures: none produced an obviously wrong result, and
two would have propagated invisibly into the final report.

---

## 1. Source

Count matrices supplied pre-processed by the course instructor, in eight folders
named `<Sex>_<Treatment>_<Replicate>`. Total 323 MB, CellRanger v3 format
(`features.tsv.gz`, not the older `genes.tsv.gz`).

Upstream: **GEO accession GSE152495**, Baker et al. 2021, *Genome Research*
31:1927–1937. Aligned by the original authors with CellRanger v3.1 against
*D. melanogaster* Release 6 (GCA_000001215.4), expected cell count 5000.

Merged dataset before QC: **88,991 cells × 17,481 genes**, identical gene count
across all eight samples.

The data is not committed to this repository — 323 MB, not ours to redistribute,
and GEO is canonical. `tools/download_data.sh` and `tools/organize_data.sh` are
retained so the dataset can be rebuilt from GEO by anyone without the
instructor's folder.

---

## 2. Correction: features.tsv.gz delimiter

**Symptom.** `sc.read_10x_mtx` failed with `KeyError: 1`.

**Diagnosis.**

```
zcat data/Female_Cocaine_1/features.tsv.gz | head -1 | awk -F'\t' '{print NF}'
-> 1        (should be 3)
```

The files are space-delimited. `read_10x_mtx` parses them with `sep="\t"` per
the CellRanger specification, so all three fields landed in column 0 and the
gene-symbol lookup `genes[1]` failed.

**Fix.** `tools/fix_features_separator.py --apply` converted all eight files to
tab-delimited, splitting on whitespace with `maxsplit=2` so the third field
("Gene Expression", which itself contains a space) stays intact. 17,481 lines
per file. Field content unchanged; only the delimiter. Originals preserved in
`data/<sample>/original/`.

---

## 3. Correction: gene symbol lost to NA parsing

**Symptom.** `adata.write()` failed with
`TypeError: Can't implicitly convert non-string objects to strings` on the `var`
index.

**Diagnosis.** Exactly one gene had a `NaN` name. Reading `features.tsv.gz` with
`keep_default_na=False` identified two candidate symbols:

| gene_id | symbol |
|---|---|
| Dmel_CG1517 | `na` (narrow abdomen) |
| Dmel_CG5842 | `nan` (nanchung) |

Pandas' default NA list contains `nan` but not lowercase `na`, so nanchung —
a TRPV channel involved in mechanosensation — became a missing value while
narrow abdomen survived as text.

**Fix.** Renamed to `CG5842` via `config.GENE_NAME_NA_FIXES`. The FlyBase symbol
is `nan`, but restoring that literal string would reintroduce the same bug at
every CSV round-trip, where a gene named `nan` is indistinguishable from a
missing value. `CG5842` is a valid identifier for the same gene and cannot be
misparsed. `na` is unaffected and retains its symbol.

**Why it mattered.** An unnamed gene never matches a marker lookup and would
appear as a blank in any results table.

---

## 4. Correction: scrambled sample labels

The most consequential of the three.

**Symptom.** `tools/verify_sample_mapping.py` reported UNCLEAR: `roX1`/`roX2`
separated the samples into two clean groups of four, but the split followed the
**treatment** field of the folder name, not the **sex** field.

**Diagnosis.** `tools/diagnose_sample_labels.py` tested two independent panels:

| Panel | Grouped by sex label | Grouped by treatment label |
|---|---|---|
| Male-biased (roX1, roX2, MSL complex) | 1.0×, overlapping | **10.6×, non-overlapping** |
| Female-biased (Sxl, Yp1–Yp3, tra) | 1.2×, overlapping | **2.1×, non-overlapping** |

Both panels tracked the treatment field, in opposite directions. The decisive
observation is the yolk proteins: `Yp1`–`Yp3` read 0.000–0.001 in all four
Sucrose-labelled samples and 0.028–0.422 in all four Cocaine-labelled ones.
Yolk proteins are transcribed only in females; no drug induces them in males.

**Cause.** The authors' published analysis code
(github.com/vshanka23/The-Drosophila-Brain-on-Cocaine-at-Single-Cell-Resolution)
lists the deposition order of the eight samples: S1–S2 female sucrose,
S3–S4 male sucrose, S5–S6 female cocaine, S7–S8 male cocaine. Pairing that order
against the supplied folder names sorted **alphabetically** reproduces the
observed marker pattern for all eight samples. The folders were built by
assigning GSM accessions in GEO order to folder names in alphabetical order —
two different orderings.

**Correction applied** (`config.SAMPLE_REMAP`, built into `scripts/01_load_data.py`):

| Folder | GEO | True sex | True treatment |
|---|---|---|---|
| Female_Cocaine_1 | S1 | Female | Sucrose |
| Female_Cocaine_2 | S2 | Female | Sucrose |
| Female_Sucrose_1 | S3 | **Male** | Sucrose |
| Female_Sucrose_2 | S4 | **Male** | Sucrose |
| Male_Cocaine_1 | S5 | **Female** | Cocaine |
| Male_Cocaine_2 | S6 | **Female** | Cocaine |
| Male_Sucrose_1 | S7 | Male | **Cocaine** |
| Male_Sucrose_2 | S8 | Male | **Cocaine** |

Four sex labels and four treatment labels changed. Original folder names are
retained in `obs['sample_folder']`.

Corrected design, balanced: Female/Cocaine 21,561 cells · Female/Sucrose 22,439 ·
Male/Cocaine 24,226 · Male/Sucrose 20,765.

### Confidence — state this distinction in the report

**Sex: confirmed.** After correction, `roX1` reads 3.898 in males vs 0.074 in
females; `Yp1`–`Yp3` are ~100× higher in females. Two independent panels,
opposite directions, non-overlapping groups.

**Treatment: inferred.** No marker gene reports whether a fly consumed cocaine.
This half rests on the deposition order in the authors' repository, not on the
data. It is supported independently: under the corrected labels the male:female
DE ratio is **8.0×** (96 vs 12 genes), and the paper reports males responding
more than females (691 vs 322, ~2.1×). That prediction played no part in
constructing the mapping.

### On the named-gene check

`tools/verify_remap.py` check 3 returned 50% agreement, but the split is
informative rather than ambiguous:

- **Up-regulated set: 5/6 agree** — IA-2 and CR34335 in both sexes, mt:lrRNA in males.
- **Down-regulated set: 0/4 agree** — roX2 and ninaE, both sexes.

If the treatment arms were swapped, *every* gene would flip, including the up
set. The up set agreeing indicates the direction is correct. Both failures have
known confounds: `roX2` is sex-linked and its variance is dominated by dosage
compensation; `ninaE` is rhodopsin, so its abundance depends on how much retinal
tissue was carried through dissection rather than on treatment. Neither is a
sound positive control.

`mt:lrRNA` at +3.034 in males is a large effect for a mitochondrial transcript.
Mitochondrial fraction is also a QC metric — check whether `pct_counts_mt`
differs by treatment before interpreting this as biology.

---

## 5. Verification record

| Check | Tool | Result |
|---|---|---|
| Delimiter converted, 8/8 files 3 columns | `fix_features_separator.py` | PASS |
| No unnamed genes survive | assertion in `01_load_data.py` | PASS |
| Design is a full 2×2×2 | `01_load_data.py` | PASS |
| Sex markers match corrected labels | `verify_remap.py` check 1 | PASS |
| Male-biased response reproduces | `verify_remap.py` check 2 | PASS, 8.0× |
| Paper's named genes | `verify_remap.py` check 3 | up set agrees, down set confounded |
| Script and notebook paths agree | both run | identical cell counts |

---

## 6. Note for the report

All three problems are silent — the pipeline runs to completion with any of them
present. The label scramble in particular would have produced a complete,
internally consistent report answering the wrong question, with the sexual
dimorphism result inverted. It was caught only because sex-specific marker genes
were checked against the supplied labels rather than assumed to agree.

Verifying inherited labels costs a minute. Not verifying them costs the project.
