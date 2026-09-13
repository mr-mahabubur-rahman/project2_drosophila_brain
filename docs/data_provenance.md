# Data provenance

## What was used in this analysis

The count matrices analysed here were **supplied pre-processed by the course
instructor**, already organized into eight sample folders:

```
data/
├── Female_Cocaine_1/   barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
├── Female_Cocaine_2/
├── Female_Sucrose_1/
├── Female_Sucrose_2/
├── Male_Cocaine_1/
├── Male_Cocaine_2/
├── Male_Sucrose_1/
└── Male_Sucrose_2/
```

Total size: 323 MB. Format: CellRanger v3 output (`features.tsv.gz`, not the
older `genes.tsv.gz`).

Upstream source: **GEO accession GSE152495**, Baker et al. 2021,
*Genome Research* 31:1927–1937. Aligned by the original authors with CellRanger
v3.1 against *D. melanogaster* Release 6 (GCA_000001215.4), expected cell count
parameter 5000.

**The data itself is not committed to this repository.** It is 323 MB, it is
not ours to redistribute, and GEO is the canonical source. `.gitignore` excludes
`data/`.

## Reproducing the data from scratch

Anyone cloning this repository without the instructor's folder can rebuild
`data/` from GEO:

```bash
bash tools/download_data.sh      # fetches GSE152495_RAW.tar
# fill in tools/sample_map.tsv from the GEO series page
bash tools/organize_data.sh      # builds the eight sample folders
```

These tools are kept in the repository deliberately, even though they were not
needed for this run. Without them, "we used GSE152495" is an assertion; with
them, it is a procedure someone else can follow.

## Label verification

Because the labels were inherited rather than assigned by us, they were verified
rather than assumed:

```bash
uv run tools/verify_sample_mapping.py
```

This checks `roX1`/`roX2`, male-specific lncRNAs of the dosage-compensation
complex. Samples labelled Male should show substantially higher expression; if
the reverse, the sex labels are inverted.

**Result:** *(record the outcome here after running it)*

Treatment labels (Cocaine vs Sucrose) have no equivalent marker and cannot be
verified computationally. They are taken as supplied. Note this as a limitation
if it becomes relevant.

## Why verify data you were given

A mislabelled sample produces a pipeline that runs flawlessly and answers the
wrong question — nothing crashes, nothing looks unusual. Inheriting labels from
a trusted source lowers the probability of error but does not eliminate it, and
the check costs one minute. Recording that you ran it is itself worth marks.
