"""
config.py — single place for every path and parameter in the pipeline.

WHY A CONFIG FILE:
During grading you will be asked "why did you choose this parameter?".
Keeping every choice in one file (instead of scattered across six scripts)
means you can point at one page and defend it. It also means a peer
reviewer can change one number and re-run without hunting through code.

PROJECT LAYOUT (config.py lives in scripts/):
    data/        the 8 CellRanger sample folders  -- gitignored
    docs/        paper.pdf, project2_guide.pdf, analysis decisions
    notebooks/   exploratory work
    peer_review/ review given to Group 2, review received, rebuttal
    reports/     draft PDF, final PDF, AI usage disclosure
    results/     checkpoints/, figures/, tables/
    scripts/     the numbered pipeline (this file)
    tools/       data download, sample mapping, environment checks
"""

from pathlib import Path

# ----------------------------------------------------------------------
# PATHS
# ----------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = PROJECT_ROOT / "data"
DOCS_DIR = PROJECT_ROOT / "docs"
NOTEBOOK_DIR = PROJECT_ROOT / "notebooks"
PEER_REVIEW_DIR = PROJECT_ROOT / "peer_review"
REPORTS_DIR = PROJECT_ROOT / "reports"
RESULTS_DIR = PROJECT_ROOT / "results"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
TOOLS_DIR = PROJECT_ROOT / "tools"

CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
FIG_DIR = RESULTS_DIR / "figures"
TABLE_DIR = RESULTS_DIR / "tables"

# Intermediate AnnData checkpoints. Each script writes one; the next reads it.
# This means a crash in step 5 does NOT force you to redo the 40-minute step 3.
# They live under results/checkpoints/ and are gitignored -- they are
# regenerable, and each is hundreds of MB.
H5AD_RAW = CHECKPOINT_DIR / "01_merged_raw.h5ad"
H5AD_QC = CHECKPOINT_DIR / "02_qc_filtered.h5ad"
H5AD_CLUSTERED = CHECKPOINT_DIR / "03_clustered.h5ad"
H5AD_ANNOTATED = CHECKPOINT_DIR / "04_annotated.h5ad"

# Reference documents (checked by tools/check_environment.py, cited in report)
PAPER_PDF = DOCS_DIR / "paper.pdf"
GUIDE_PDF = DOCS_DIR / "project2_guide.pdf"

for _d in (RESULTS_DIR, CHECKPOINT_DIR, FIG_DIR, TABLE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------
# SAMPLE LABEL CORRECTION  —  read docs/data_provenance.md before changing
# ----------------------------------------------------------------------
# The supplied folder names DO NOT match their contents.
#
# Evidence: sex-specific markers separate the eight samples cleanly into two
# groups of four, but the split follows the TREATMENT field of the folder name,
# not the SEX field.
#     roX1/roX2 (male-specific)  high in the four *_Sucrose_* folders
#     Yp1-Yp3, Sxl (female)      high in the four *_Cocaine_* folders
# Yolk proteins are transcribed only in females and read 0.000-0.001 in all four
# Sucrose-labelled samples, so this cannot be a drug effect.
#
# Cause: the authors' published R code
# (github.com/vshanka23/The-Drosophila-Brain-on-Cocaine-at-Single-Cell-Resolution)
# lists the deposition order of the eight samples:
#     S1 Female_sucrose_R1   S5 Female_cocaine_R1
#     S2 Female_sucrose_R2   S6 Female_cocaine_R2
#     S3 Male_sucrose_R1     S7 Male_cocaine_R1
#     S4 Male_sucrose_R2     S8 Male_cocaine_R2
# Pairing that order against the folder names sorted ALPHABETICALLY reproduces
# the observed marker pattern for all eight samples. The folders were built by
# assigning GSM accessions in GEO order to folder names in alphabetical order.
#
# CONFIDENCE — state this distinction in your Methods:
#   sex       CONFIRMED. Two independent marker panels, opposite directions,
#             non-overlapping groups.
#   treatment INFERRED from the deposition order. No marker reports whether a
#             fly ate cocaine. Supported independently: under these labels the
#             male:female DE ratio is ~8x, and the paper reports males
#             responding more than females (691 vs 322, ~2.1x).
#
# folder name -> (sex, treatment, replicate, GEO position)
SAMPLE_REMAP = {
    "Female_Cocaine_1": ("Female", "Sucrose", "1", "S1"),
    "Female_Cocaine_2": ("Female", "Sucrose", "2", "S2"),
    "Female_Sucrose_1": ("Male",   "Sucrose", "1", "S3"),
    "Female_Sucrose_2": ("Male",   "Sucrose", "2", "S4"),
    "Male_Cocaine_1":   ("Female", "Cocaine", "1", "S5"),
    "Male_Cocaine_2":   ("Female", "Cocaine", "2", "S6"),
    "Male_Sucrose_1":   ("Male",   "Cocaine", "1", "S7"),
    "Male_Sucrose_2":   ("Male",   "Cocaine", "2", "S8"),
}

# Set False only to reproduce the original (incorrect) labelling for comparison.
APPLY_SAMPLE_REMAP = True

# Gene symbols pandas misreads as missing values when parsing features.tsv.
# 'nan' is nanchung (Dmel_CG5842), a TRPV channel. Left unhandled it becomes
# NaN in var_names, which breaks h5ad writing and every marker lookup.
# 'na' (narrow abdomen, Dmel_CG1517) is NOT affected -- lowercase 'na' is not
# in pandas' default NA list.
GENE_NAME_NA_FIXES = {
    "Dmel_CG5842": "CG5842",   # 'nan' / nanchung
}

# ----------------------------------------------------------------------
# QUALITY CONTROL
# ----------------------------------------------------------------------
# The paper (Baker et al. 2021, Methods) filtered cells to 300-2500 detected
# genes using Seurat's recommendation, and dropped genes seen in <5 cells.
# The project guide relaxes this to min_genes=200 / min_cells=3.
# We follow the PAPER, because Research Question 1 asks us to reproduce their
# 36 clusters -- using their thresholds removes one source of divergence.
# Set USE_PAPER_THRESHOLDS = False to fall back to the guide's numbers.
USE_PAPER_THRESHOLDS = True

MIN_GENES_PER_CELL = 300 if USE_PAPER_THRESHOLDS else 200
MAX_GENES_PER_CELL = 2500 if USE_PAPER_THRESHOLDS else None  # None = no upper cut
MIN_CELLS_PER_GENE = 5 if USE_PAPER_THRESHOLDS else 3

# Mitochondrial cut-off. The paper does NOT apply one; the guide asks for <10%.
# We keep it, because high mt% marks stressed/dying cells from the dissociation
# protocol, and brain dissociation is harsh. Inspect the QC violins before
# committing to this number -- if your distribution peaks above 10 you are
# throwing away real cells.
MAX_PCT_MT = 10.0

# ----------------------------------------------------------------------
# NORMALIZATION / FEATURE SELECTION
# ----------------------------------------------------------------------
TARGET_SUM = 1e4          # counts-per-10k, the scanpy convention
N_TOP_GENES = 2000        # HVGs, per the guide
HVG_BATCH_KEY = "sample"  # select HVGs per-sample then rank by consistency,
                          # so genes driven by one bad sample don't dominate

# ---- MEMORY PROFILE: tuned for a 16 GB machine ----------------------
# On WSL2 this is the number that matters, and it is NOT your host RAM.
# WSL2 defaults to giving Linux half the host RAM (so 8 GB of your 16), which
# is not enough for step 03. Fix it with a .wslconfig file -- see README §1.
# Check what Linux can actually see:   free -g
#
# regress_out fits a linear model per gene across all cells. We subset to HVGs
# BEFORE regressing, so it operates on 2,000 genes rather than ~13,000. At
# ~80k cells that is a dense 80,000 x 2,000 float32 array (~640 MB) -- fits in
# 16 GB, but expect 15-30 minutes.
#
# If step 03 gets killed by the OOM killer, turn this off FIRST. It is the
# largest single saving, and modern scanpy tutorials have dropped the step
# entirely -- the biological conclusions barely move.
REGRESS_OUT = True
REGRESS_VARS = ["total_counts", "pct_counts_mt"]
SCALE_MAX_VALUE = 10

# Last resort, only after REGRESS_OUT=False has failed: cap the cell count
# before step 03. None = use all cells. 40000 roughly halves the footprint.
# If you use this, you MUST report it in Methods -- it changes cluster counts
# and DE power, and a reviewer comparing your cell count to the paper's 86,224
# will ask.
SUBSAMPLE_N_CELLS = None

# ----------------------------------------------------------------------
# DIMENSIONALITY REDUCTION & CLUSTERING
# ----------------------------------------------------------------------
N_PCS = 30            # inspect the elbow in the PCA variance plot before trusting this
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 0.8   # the paper's value -- 36 clusters
LEIDEN_KEY = "leiden_res_0.8"
RANDOM_SEED = 0

# Research Question 1 asks whether YOU reproduce 36 clusters. Scanpy+Leiden
# will not land on exactly 36 at res 0.8 -- the paper used Seurat's SNN +
# SCTransform, a different algorithm on differently normalized data.
# This sweep lets you report the resolution that DOES give ~36, which is a
# far better answer than "we got 29 and moved on".
RESOLUTION_SWEEP = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
TARGET_N_CLUSTERS = 36

# Optional batch integration across the 8 samples (needs `harmonypy`).
# The paper reports no single cluster was dominated by one sample, so heavy
# correction is probably unnecessary -- but check the UMAP-by-sample first.
USE_HARMONY = False
HARMONY_BATCH_KEY = "sample"

# ----------------------------------------------------------------------
# DIFFERENTIAL EXPRESSION
# ----------------------------------------------------------------------
DE_METHOD = "wilcoxon"
# NOTE ON UNITS: scanpy reports log2 fold change. The paper reports NATURAL log
# (|log_e FC| > 1). So the paper's threshold is equivalent to |log2FC| > 1.44.
# The guide asks for |log2FC| > 1.0, which is the LOOSER cut. We report both.
LOG2FC_THRESHOLD = 1.0
PAPER_EQUIVALENT_LOG2FC = 1.4427  # = 1 / ln(2)
PADJ_THRESHOLD = 0.05

# Per-cluster DE is skipped for clusters with fewer than this many cells in
# either treatment arm -- Wilcoxon on 4 cells produces noise, not biology.
MIN_CELLS_PER_GROUP_FOR_DE = 20

# ----------------------------------------------------------------------
# PATHWAY ENRICHMENT
# ----------------------------------------------------------------------
# Enrichr keeps SEPARATE library names per organism. "KEGG_2021_Human" (in the
# project guide) does NOT exist in the Fly modality and will silently return
# nothing useful. 06_pathway_enrichment.py auto-detects what is actually
# available and falls back to these.
ENRICHR_ORGANISM = "fly"
ENRICHR_LIBRARIES_PREFERRED = [
    "GO_Biological_Process_2018",
    "GO_Molecular_Function_2018",
    "KEGG_2019",
]
N_GENES_FOR_ENRICHMENT = 150

# ----------------------------------------------------------------------
# MARKER GENES — RESOLVED against this dataset's reference annotation
# ----------------------------------------------------------------------
# This reference substitutes vertebrate ortholog names for five fly symbols.
# Verified by gene ID with tools/resolve_markers.py; worksheet in
# results/tables/marker_resolution.csv. The symbols below are the ones that
# actually appear in adata.var_names -- the true fly symbol is in the comment.
#
# CASE SENSITIVITY WARNING: in this reference `trh` (lowercase) is
# Dmel_CG42865 = trachealess, a tracheal transcription factor. The serotonin
# synthesis enzyme Trh (CG9122) is called `Trhn` here. A case-insensitive
# marker match finds trachealess, reports success, and would label a cluster
# "serotonergic" on the strength of a tracheal gene. Never match markers
# case-insensitively in Drosophila.
CANONICAL_MARKERS = {
    "Glia (pan)": ["repo"],                                  # CG31240
    "Astrocyte-like glia": ["alrm", "Eaat1", "Gat"],          # CG11910, CG3747, CG1732
    "Surface glia / BBB": ["moody", "Mdr65", "Tret1"],        # CG4322, CG10181, CG30035 (= Tret1-1)
    "Cortex/ensheathing glia": ["wrapper", "zyd"],            # CG10382, CG2893
    "Neuron (pan)": ["elav", "nSyb", "brp", "Syt1"],          # CG4262, CG17248, CG42344, CG3139
    "Kenyon cells (MB)": ["ey", "Fas2", "Adcy1", "Pde4",      # CG1464, CG3665, CG9533 (= rut), CG32498 (= dnc)
                          "sNPF", "Dop1R2"],                  # CG13968, CG18741
    "Cholinergic": ["VAChT", "ChAT"],                         # CG32848, CG12345
    "GABAergic": ["Gad1", "VGAT"],                            # CG14994, CG8394
    "Glutamatergic": ["VGlut1"],                              # CG9887 (= VGlut; NOT the vertebrate paralog)
    "Dopaminergic": ["ple", "DAT", "Ddc"],                    # CG10118, CG8380, CG10697
    "Serotonergic": ["SerT", "Trhn"],                         # CG4545, CG9122 (= Trh)
    "Octopaminergic": ["Tdc2", "Tbh"],                        # CG30446, CG1543
    "Photoreceptor": ["ninaE", "Rh3", "Rh4", "trp"],          # CG4550, CG10888, CG9668, CG7875
    "Monoaminergic (vesicular)": ["Vmat"],                    # CG33528
}

# Reference symbol -> true fly symbol. Use these when labelling figures and
# writing the report, so a reader is not left wondering why a Drosophila paper
# cites VGlut1. Note this reference also contains a separate VGlut2 (CG4288) --
# these are NOT the vertebrate paralogs.
REFERENCE_SYMBOL_ALIASES = {
    "VGlut1": "VGlut",
    "Adcy1": "rut",
    "Pde4": "dnc",
    "Trhn": "Trh",
    "Tret1": "Tret1-1",
}

# Sex-linked controls. roX1/roX2 are male-specific lncRNAs -- used by
# tools/verify_sample_mapping.py to catch an inverted sex label.
SEX_CONTROL_GENES = ["roX1", "roX2", "lncRNA:roX1", "lncRNA:roX2"]

# Genes the paper highlights as globally cocaine-responsive (Results, TopKLists).
# Use these as a positive control: if none of them move in your DE, something
# upstream is wrong.
PAPER_CORE_RESPONSE_GENES = {
    "up": ["Rpl41", "IA-2", "CR34335", "mt:lrRNA"],
    "down": ["roX2", "ninaE"],
}

# Genes the paper singles out in the Discussion -- worth checking by name.
PAPER_DISCUSSION_GENES = [
    "rut", "cpx", "Ddc", "Dop2R", "slo", "Rgk1", "jdp",
    "ogre", "Nrg", "moody", "Eaat1", "CG3168",
    "GstE12", "GstE14", "se", "Ssadh", "SPARC", "bdl", "Tsp", "Fas2",
    "Aldh", "GluRIA", "GluRIB", "Vmat",
]

# Clusters the paper reports as the strongest responders. Their "C11" numbering
# is NOT your Leiden numbering -- match by marker genes, never by number.
PAPER_TOP_RESPONDER_CLUSTERS = {
    "C11": "Kenyon cells type 1",
    "C16": "unannotated (antennal/optic lobe mix)",
    "C17": "astrocytes",
    "C22": "surface glia and fat body",
}
