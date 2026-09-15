# Repository description and README

---

## Short description (GitHub "About" field — 350 char limit)

Pick one:

**Option A — factual**
> Independent Scanpy reanalysis of GSE152495 (Baker et al. 2021, *Drosophila* brain single-cell cocaine response). Reproduces the published cell count to 0.05% and tests the differential expression results against five statistical controls, including a permutation comparison against alternative sample groupings.

**Option B — leads with the finding**
> Reanalysis of a published *Drosophila* single-cell cocaine dataset (GSE152495). The design-level findings reproduce; the reported female response falls below its own background when compared against alternative groupings of the same samples.

**Option C — short**
> Independent reanalysis and statistical verification of GSE152495 — single-cell transcriptomics of the cocaine-exposed *Drosophila* brain.

**Topics / tags:**
`single-cell-rna-seq` · `scanpy` · `drosophila` · `reproducibility` ·
`bioinformatics` · `differential-expression` · `reanalysis` · `geo-dataset`

---

## README.md

```markdown
# Reanalysis of the cocaine-exposed *Drosophila* brain (GSE152495)

Independent Scanpy-based reanalysis of the single-cell transcriptomic atlas
published by **Baker et al. (2021)**, *Genome Research* 31:1927–1937
([doi:10.1101/gr.268037.120](https://doi.org/10.1101/gr.268037.120)).

The aim was twofold: to test whether the published findings **reproduce** under a
different pipeline, and to test whether they **survive** statistical controls that
the original analysis did not report.

---

## Summary of findings

| Published finding | This reanalysis |
|---|---|
| 86,224 cells after QC | 86,177 — 0.05% difference |
| 36 clusters at resolution 0.8 | 30; reaches 36 at ~1.33 |
| Cluster count plateaus at 0.8 | Not reproduced — monotonic increase |
| Kenyon cells among top responders | 1st on raw count, 7th at matched power |
| Male-biased response, 2.15× | 0.97× to 6.43× depending on the metric |
| Female response, 322 genes | Falls below its own background |

The central result is a permutation control. Within each sex, four samples permit
two null contrasts besides the true treatment split:

| | Treatment (real) | Replicate axis (null) | Diagonal (null) |
|---|---|---|---|
| Male | **90** | 7 | 14 |
| Female | **14** | 39 | 35 |

In males the treatment effect exceeds both nulls by 6–13×. In females both nulls
exceed it, and 64% of the 14 female genes also appear in a null contrast. The
design-level conclusions of the original study reproduce; the female response is
not resolvable above between-sample variation at this replicate number.

---

## Three data problems found before analysis

Each would have propagated silently:

1. **`features.tsv.gz` was space- rather than tab-delimited**, so
   `scanpy.read_10x_mtx` could not assign gene symbols.
2. **The gene *nanchung* (`nan`, Dmel_CG5842) is parsed by pandas as a missing
   value**, leaving one gene unnamed and blocking HDF5 serialisation.
3. **Sample folder labels did not match their contents.** Sex-specific markers
   separated the eight samples by the *treatment* field of the folder name, not
   the *sex* field. Identities were reconstructed from the deposition order in
   the authors' published code and verified against two independent marker
   panels.

The reference also substitutes vertebrate ortholog names for five fly symbols
(`VGlut1`=VGlut, `Adcy1`=*rut*, `Pde4`=*dnc*, `Trhn`=Trh, `Tret1`=Tret1-1).
Note that `trh` lowercase in this reference is *trachealess*, an unrelated gene —
marker matching must be case-sensitive.

---

## Repository layout

```
├── data/          8 CellRanger sample folders          [gitignored]
├── docs/          reference PDFs, provenance, decisions log
├── notebooks/     the pipeline as JupyterLab notebooks
├── peer_review/   review checklist and rebuttal templates
├── reports/       manuscript drafts and AI usage disclosure
├── results/
│   ├── checkpoints/   .h5ad between stages              [gitignored]
│   ├── figures/       all report figures
│   └── tables/        CSV outputs
├── scripts/       config.py + numbered pipeline (01–10)
└── tools/         data download, label verification, figure generation
```

`scripts/config.py` holds every analysis parameter with the reasoning for each,
so a reviewer can change one number and re-run rather than searching the code.

---

## Reproducing this analysis

```bash
# environment
uv venv --python 3.11 && source .venv/bin/activate
uv pip install -r requirements.txt
python tools/check_environment.py

# data — matrices are not redistributed here
bash tools/download_data.sh          # fetches GSE152495 from GEO
# fill in tools/sample_map.tsv, then
bash tools/organize_data.sh

# pipeline
python scripts/01_load_data.py
python tools/verify_sample_mapping.py    # confirms labels using roX1/roX2
python scripts/02_qc_filter.py
python scripts/03_normalize_cluster.py
python scripts/04_annotate_clusters.py   # fill the worksheet, then re-run
python scripts/05_de_analysis.py

# controls
python scripts/05b_pseudobulk_check.py   # replicate-level validation
python scripts/06_pathway_enrichment.py
python scripts/07_interaction_test.py    # sex × treatment interaction
python scripts/08_sensitivity.py         # gene-class exclusion, size matching
python scripts/09_reproduce_fig3.py      # Baker et al. Figure 3A–C
python scripts/10_permutation_control.py # DE against alternative groupings
```

Every stage writes an `.h5ad` checkpoint, so a failure at step 5 does not
require redoing step 3. All stochastic steps use a fixed seed.

---

## Statistical controls

| Script | Question | Result |
|---|---|---|
| `05b` | Do cell-level results hold at replicate level? | 100% direction agreement, ρ = 0.78 |
| `07` | Does any gene respond *differently* by sex? | 0 genes at FDR < 0.05 (4 residual df) |
| `08A` | Do mitochondrial genes drive the male bias? | No — pooled ratio unchanged |
| `08B` | Is the response ranking power or biology? | Kenyon cells 1st → 7th at matched power |
| `10` | Is the effect larger than its own background? | Male yes; female no |

---

## Environment

Python 3.11.15 · scanpy 1.11.5 · anndata 0.12.19 · numpy 2.4.6 · pandas 2.3.3 ·
scipy 1.17.1 · scikit-learn 1.9.1 · leidenalg 0.12.0 · umap-learn 0.5.12 ·
gseapy 1.3.1 · statsmodels. Exact versions in `requirements.lock.txt`.

---

## Data availability

Raw data: **GEO accession [GSE152495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152495)**.

Count matrices are not redistributed in this repository. `tools/download_data.sh`
and `tools/organize_data.sh` rebuild `data/` from GEO.

---

## Limitations

Stated fully in the report; the main ones:

- **Treatment labels are inferred**, not verified. Sex was confirmed by marker
  genes; treatment assignment rests on the deposition order in the authors'
  published code. No computational check can verify which flies consumed cocaine.
- **Two null contrasts per sex** is all four samples allow. This gives an
  order-of-magnitude comparison, not a *p*-value.
- **Absence of evidence is not evidence of absence.** That the female response is
  not resolvable above background does not show that female flies do not respond.
- **Six clusters (15.7% of cells) are unannotated**, including the one ranking
  highest under size-matched analysis.
- **One cluster was excluded post hoc** after its depth imbalance was observed.

---

## Citation

If you use this code, please cite the original study:

> Baker BM, Mokashi SS, Shankar V, Hatfield JS, Hannah RC, Mackay TFC, Anholt RRH.
> 2021. The *Drosophila* brain on cocaine at single-cell resolution.
> *Genome Research* 31: 1927–1937. doi:10.1101/gr.268037.120

---

## Acknowledgement

Analysis performed as a self-directed reanalysis project. AI assistance was used
during pipeline development and is documented in
`reports/AI_USAGE_DISCLOSURE.md`, including prompts and validation steps.
```

---

## Notes before you publish

- **Check the `git status --short` output** before pushing to confirm `data/` and
  `results/checkpoints/` are excluded. You should be pushing roughly 50–60 files,
  not thousands.
- **The AI acknowledgement at the end is optional for a public repo** but
  consistent with your course policy, and stating it plainly is better than
  having someone infer it from the commit history.
- **Consider adding a LICENSE.** MIT is typical for analysis code. The data is
  not yours to license, which is another reason not to commit it.
- **If you make the repo public**, the README claims are checkable against your
  own tables — make sure the numbers quoted here match `results/tables/`.
