# Reanalysis of the cocaine-exposed *Drosophila* brain (GSE152495)

Independent Scanpy-based reanalysis of the single-cell transcriptomic atlas
published by **Baker et al. (2021)**, *Genome Research* 31:1927-1937
([doi:10.1101/gr.268037.120](https://doi.org/10.1101/gr.268037.120)).

Two aims: to test whether the published findings **reproduce** under a different
pipeline, and whether they **survive** statistical controls the original analysis
did not report.

![Graphical abstract](results/figures/fig0_graphical_abstract.png)

## Summary of findings

| Published finding | This reanalysis |
|---|---|
| 86,224 cells after QC | **86,177** - 0.05% difference |
| 36 clusters at resolution 0.8 | **30**; reaches 36 at ~1.33 |
| Cluster count plateaus at 0.8 | **Not reproduced** - monotonic increase |
| Kenyon cells among top responders | 1st raw, **7th at matched power** |
| Male-biased response, 2.15x | **0.97x to 6.43x** depending on the metric |
| Female response, 322 genes | **Falls below its own background** |

### The central result

Within each sex, four samples permit two null contrasts besides the true
treatment split, using identical cells, test and thresholds:

| | Treatment (real) | Replicate axis (null) | Diagonal (null) |
|---|---|---|---|
| **Male** | **90** | 7 | 14 |
| **Female** | **14** | 39 | 35 |

In males the treatment effect exceeds both nulls by 6-13x. In females both nulls
exceed it, and 64% of the 14 female genes also appear in a null contrast -
including the yolk proteins Yp1-Yp3, independently identified as ambient RNA.

**The design-level conclusions reproduce. The female response is not resolvable
above between-sample variation at this replicate number.**

## Three data problems found before analysis

1. **features.tsv.gz was space- rather than tab-delimited**, so scanpy could not
   assign gene symbols.
2. **The gene nanchung (`nan`, Dmel_CG5842) is parsed by pandas as a missing
   value**, leaving one gene unnamed.
3. **Sample folder labels did not match their contents.** Sex markers separated
   the samples by the treatment field, not the sex field. Identities were
   reconstructed from the deposition order in the authors' published code.

The reference also substitutes vertebrate ortholog names for five fly symbols.
Note `trh` lowercase here is *trachealess*, an unrelated gene - marker matching
must be case-sensitive.

## Statistical controls

| Script | Question | Result |
|---|---|---|
| 05b | Do cell-level results hold at replicate level? | 100% direction agreement |
| 07 | Does any gene respond differently by sex? | 0 genes at FDR < 0.05 |
| 08A | Do mitochondrial genes drive the male bias? | No - ratio unchanged |
| 08B | Is the ranking power or biology? | Kenyon cells 1st to 7th |
| 10 | Is the effect larger than its own background? | Male yes; female no |

## Repository layout

    data/          8 CellRanger sample folders       [not tracked]
    docs/          reference PDFs, provenance, setup guide
    notebooks/     the pipeline as JupyterLab notebooks
    peer_review/   review checklist and rebuttal templates
    reports/       final report (docx + pdf)
    results/       figures, tables, checkpoints
    scripts/       config.py + numbered pipeline (01-10)
    tools/         data download, label verification, figures

Setup instructions are in [docs/SETUP.md](docs/SETUP.md).

## Data availability

Raw data: GEO accession
[GSE152495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152495).
Count matrices are not redistributed here.

## Limitations

- **Treatment labels are inferred, not verified.** Sex was confirmed by markers.
- **Two null contrasts per sex** is all four samples allow.
- **Absence of evidence is not evidence of absence** for the female response.
- **Six clusters (15.7% of cells) are unannotated.**
- **One cluster was excluded post hoc** after its depth imbalance was observed.

## Citation

> Baker BM, Mokashi SS, Shankar V, Hatfield JS, Hannah RC, Mackay TFC, Anholt
> RRH. 2021. The Drosophila brain on cocaine at single-cell resolution.
> Genome Research 31: 1927-1937. doi:10.1101/gr.268037.120

## Acknowledgement

Self-directed reanalysis project. AI assistance is documented in
reports/AI_USAGE_DISCLOSURE.md.
