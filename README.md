# Project 2 — The *Drosophila* Brain on Cocaine (scRNA-seq reanalysis)

Reanalysis of **GSE152495** (Baker et al. 2021, *Genome Research* 31:1927–1937)
using Scanpy. Eight 10x samples: Sex (Female/Male) × Treatment (Cocaine/Sucrose)
× 2 replicates.

---

## Layout

```
project2_drosophila_brain/
├── data/          8 CellRanger sample folders          [gitignored]
├── docs/          paper.pdf, project2_guide.pdf, analysis_decisions.md
├── notebooks/     the same pipeline as JupyterLab notebooks
├── peer_review/   review given, review received, rebuttal
├── reports/       draft PDF, final PDF, AI usage disclosure, report outline
├── results/
│   ├── checkpoints/   .h5ad between pipeline stages    [gitignored]
│   ├── figures/       report figures
│   └── tables/        CSV outputs
├── scripts/       config.py + the numbered pipeline
└── tools/         data download, sample mapping, environment checks
```

The split that matters: **`scripts/` is the pipeline a reviewer reruns;
`tools/` is one-time setup; `notebooks/` is for looking, not deciding.**
Analysis decisions belong in `scripts/config.py`, so the run is reproducible
from the command line rather than by executing cells in the right order.

---

## 1. Environment

### 1a. Give WSL enough memory — do this first

On a 16 GB PC, **WSL2 hands Linux only 8 GB by default** (half the host, capped
at 8). That is not enough for step 03, and the failure mode is the kernel being
killed with no useful message.

In Windows, create `C:\Users\<you>\.wslconfig`:

```ini
[wsl2]
memory=12GB
swap=8GB
```

Then in PowerShell:

```powershell
wsl --shutdown
```

Reopen your terminal and confirm Linux can actually see it:

```bash
free -g          # "total" should now read ~12, not ~7
```

Leaving 4 GB for Windows is deliberate — starving the host makes everything,
including WSL, slower.

### 1b. Install

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.bashrc

cd project2_drosophila_brain
uv venv --python 3.11
source .venv/bin/activate
uv pip install -r requirements.txt
uv pip freeze > requirements.lock.txt     # commit this
```

### 1c. Verify before running anything long

```bash
uv run tools/check_environment.py
```

Checks the Leiden backend (`igraph` + `leidenalg`, which scanpy does not install
itself and which fails unhelpfully at step 03), the folder layout, the data, and
the RAM Linux can actually see.

**Work inside the WSL filesystem** (`/home/<you>/...`), not `/mnt/c/...`.
Reading thousands of matrix entries across the Windows boundary is roughly 10×
slower, and `check_environment.py` will warn you if you forget.

---

## 2. Data

```bash
bash tools/download_data.sh      # GSE152495_RAW.tar from GEO, unpacked
# fill in tools/sample_map.tsv by hand — see that file's header
bash tools/organize_data.sh      # builds data/<Sample_Name>/
```

Then, after step 01:

```bash
uv run tools/verify_sample_mapping.py
```

This is not optional. A swapped sex or treatment label produces a pipeline that
runs flawlessly and answers the wrong question — nothing crashes, nothing looks
odd. The check uses `roX1`/`roX2`, male-specific lncRNAs, to catch an inverted
sex mapping. Treatment labels have no such marker and must be verified by hand
against the GEO page and Supplemental Table S2.

---

## 3. Run order

| Script | Does | Roughly |
|---|---|---|
| `scripts/01_load_data.py` | Load 8 samples, attach metadata, merge | 2–5 min |
| `tools/verify_sample_mapping.py` | Confirm labels are not inverted | 1 min |
| `scripts/02_qc_filter.py` | QC metrics, plots, filtering → **Fig 1** | 3–8 min |
| `scripts/03_normalize_cluster.py` | Normalize, HVG, PCA, UMAP, Leiden → **Fig 2** | 20–60 min |
| `scripts/04_annotate_clusters.py` | Markers + cell types → **Fig 3** | 10–20 min |
| `scripts/05_de_analysis.py` | DE: global, per-cluster, sex overlap → **Fig 4** | 15–40 min |
| `scripts/05b_pseudobulk_check.py` | Replicate-level validation (recommended) | 5 min |
| `scripts/06_pathway_enrichment.py` | GSEApy/Enrichr → **Fig 5** | 5–15 min |

```bash
uv run scripts/01_load_data.py
uv run tools/verify_sample_mapping.py
uv run scripts/02_qc_filter.py
uv run scripts/03_normalize_cluster.py
uv run scripts/04_annotate_clusters.py
# --- fill results/tables/annotation_worksheet.csv,
#     save as annotation_worksheet_FILLED.csv, re-run 04 ---
uv run scripts/04_annotate_clusters.py
uv run scripts/05_de_analysis.py
uv run scripts/05b_pseudobulk_check.py
uv run scripts/06_pathway_enrichment.py
```

Each script writes a checkpoint to `results/checkpoints/`, so a crash in step 5
never forces you to redo step 3.

---

## 3b. The JupyterLab path

Every stage also exists as a notebook in `notebooks/`, covering the same
operations with inspection points between them.

```bash
python -m ipykernel install --user --name project2 --display-name "Python 3 (project2)"
uv run jupyter lab
```

Then work through `01_load_data.ipynb` → `02_qc_filter.ipynb` → … in order.
Select the **Python 3 (project2)** kernel, not the system Python.

### How the two paths stay in sync

Both import every parameter from `scripts/config.py`. The notebooks never
hard-code a threshold or a resolution. So the two paths cannot disagree about an
analysis decision — change it in `config.py` and both follow.

The notebooks are generated by `tools/build_notebooks.py`. Edit that file and
re-run it rather than hand-editing the `.ipynb` JSON, whose git diffs are
unreadable and effectively un-reviewable.

### Which path for what

| | Use |
|---|---|
| `scripts/` | The reproducible run. One command per stage, no hidden state, no execution order to get wrong. **This is what you tag for the milestones and what a peer reviewer runs.** |
| `notebooks/` | Learning the operations, inspecting QC distributions, working through cluster annotation gene by gene, and debugging. |

A notebook can be run out of order or with a stale variable still in memory, and
neither leaves a trace — which is exactly why the scripts are the submitted
artefact. Use the notebooks to *look*; use `config.py` to *decide*.

### Memory note for notebooks (16 GB)

**Restart the kernel between notebooks.** A kernel still holding `adata` from
notebook 02 is the usual cause of an OOM kill in notebook 03. Close other
notebooks, and prefer `Restart Kernel and Run All Cells` over re-running cells
piecemeal.

If step 03 is killed anyway, set `REGRESS_OUT = False` in `config.py` first —
it is the largest single saving. `SUBSAMPLE_N_CELLS` is the last resort, and
using it must be declared in your Methods.

---

## 4. Where you have to think, not just run

Three places the pipeline deliberately stops and hands the problem back to you:

**Cluster annotation (step 4).** The script writes a worksheet with top markers
per cluster and leaves the cell-type column blank. Filling it in is the actual
scientific work. The paper's C11 is Kenyon cells; **your** cluster 11 is almost
certainly something else, because Leiden numbers clusters by size and your
clustering is not theirs. Match on marker genes, never on number.

**QC thresholds (step 2).** The script plots before it cuts. Look at the violins
before accepting the defaults in `config.py`.

**Number of PCs (step 3).** `N_PCS=30` is a starting point. Check the elbow in
`results/figures/pca_variance_03_pca_variance.png`.

---

## 5. Things that will differ from the paper, and why

You will not reproduce their numbers exactly. Saying so clearly is worth more
marks than pretending otherwise.

| They did | We do | Consequence |
|---|---|---|
| Seurat v3, SCTransform | Scanpy, CPM + log1p | Different variance structure → different cluster boundaries |
| SNN clustering | Leiden | Related but not identical graph partitioning |
| MAST for DE | Wilcoxon | Different power, different gene counts |
| `\|log_e FC\| > 1` | `\|log2FC\| > 1` | Their cut is *stricter* (= `\|log2FC\| > 1.44`). Step 05 reports both. |
| Bonferroni | Benjamini–Hochberg | BH is less conservative → more hits |

Expect roughly 30–40 clusters rather than exactly 36, and different absolute DE
counts. What **should** reproduce is the pattern: a strong male bias in DE
burden, with glia and Kenyon cells among the top responders. Report the ratio
and the ranking, not the raw counts.

---

## 6. Milestones

```bash
git tag -a v1.0-peerreview -m "Milestone 1: draft report"
git push origin v1.0-peerreview

git tag -a v2.0-final -m "Milestone 3: final submission"
git push origin v2.0-final
```

Commit code, `requirements.lock.txt`, `results/tables/`, `results/figures/`,
`docs/analysis_decisions.md`, and the reports. Never commit `data/` or the
`.h5ad` checkpoints — `.gitignore` handles both.

See `peer_review/README.md` for the Milestone 2 review checklist and the
Milestone 3 rebuttal format, and `reports/report_outline.md` for the mapping
from pipeline outputs to report sections.
