# docs/

Reference material and the decisions log.

| File | What it is |
|---|---|
| `paper.pdf` | Baker et al. 2021, *Genome Research* 31:1927–1937 — the primary reference |
| `project2_guide.pdf` | Course project guide: timeline, research questions, deliverables |
| `analysis_decisions.md` | Running log of every choice and its justification |

---

# analysis_decisions.md — starter

Keep this current as you work. At grading you must justify every parameter, and
a log written at the time is more accurate than one reconstructed afterwards.

## Decisions carried over from setup

| Decision | Chosen | Why | Guide said |
|---|---|---|---|
| QC gene range | 300–2500 per cell | Paper's Methods; matching their filter removes one source of divergence from their 36 clusters | 200, no upper bound |
| Min cells per gene | 5 | Paper's Methods | 3 |
| Mitochondrial cut | < 10% | Guide's value; brain dissociation is harsh and high mt% marks dying cells. Paper applied none | < 10% |
| Upper gene bound | 2500 | Crude doublet filter — droplets with more detected genes than a single fly neuron plausibly expresses | not specified |
| Normalization | CPM 1e4 + log1p | Scanpy standard. Paper used SCTransform — a genuinely different variance model, and part of why cluster counts differ | same |
| HVGs | 2000, `batch_key='sample'` | Per-sample ranking stops one anomalous sample driving selection | 2000, `batch_key='sample'` |
| HVG subset before scaling | Yes | Makes `regress_out` tractable at ~80k cells (2,000 genes not ~13,000). `adata.raw` keeps all genes for DE | not specified |
| Leiden resolution | 0.8, plus a sweep 0.4–1.4 | 0.8 is the paper's value; the sweep reproduces their stability argument and handles the near-certain mismatch honestly | 0.8 |
| DE test | Wilcoxon | Scanpy default; non-parametric, appropriate for zero-inflated counts. Paper used MAST | Wilcoxon |
| DE threshold | \|log2FC\| > 1, BH p < 0.05; **paper equivalent 1.44 also reported** | Scanpy reports log2, the paper natural log. Reporting both prevents a false comparison | \|log2FC\| > 1 |
| Enrichr libraries | Auto-detected from the Fly catalogue | The guide's `KEGG_2021_Human` does not exist in the Fly modality and returns nothing useful | GO_BP_2023, KEGG_2021_Human |
**Cluster 11 (res 0.8), n=3,032.** Neuronal (elav/nSyb high, repo low), top
markers Imp / futsch / mbl / mamo / CG31345 — corresponds to the paper's
mbl/Imp/CG31345 cluster. 70.8% female overall; 74.5% female among depth-matched
cells (1200–1600 genes) against a 56.4% background, so the sex bias is
biological, not a capture artefact. Fat body contamination ruled out (Yp1–Yp3
at background). Note: median depth 1,376 genes vs 835 elsewhere, and several
respiratory-chain and ribosomal genes appear among its markers — the cluster
boundary is partly depth-influenced.

## Decisions still to make

- [ ] `N_PCS` — currently 30. Check the elbow in `pca_variance_03_pca_variance.png`.
- [ ] `REGRESS_OUT` — keep True or drop, depending on available RAM.
- [ ] `USE_HARMONY` — only if the UMAP-by-sample shows sample-driven clusters.
- [ ] Cluster annotations — the worksheet in `results/tables/`.
- [ ] Whether to report the resolution-0.8 clustering or the one nearest 36 as primary.

## Log

*(Append dated entries as you go: what you changed, what you observed, what you concluded.)*

**2026-09-12** — Project scaffolded. Pipeline written, data not yet downloaded.
