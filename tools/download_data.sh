#!/usr/bin/env bash
# download_data.sh — fetch GSE152495 from GEO and organize into data/.
#
# The pipeline expects:
#   data/Female_Cocaine_1/{barcodes.tsv.gz,features.tsv.gz,matrix.mtx.gz}
#   ... eight folders, one per sample.
#
# GEO does not hand you that layout. Supplementary files arrive in one tar
# with GSM-prefixed names, and the exact naming varies between submissions.
# So: download, unpack, INSPECT, then map GSM IDs to readable names.
#
# Usage:  bash tools/download_data.sh
# Size:   roughly 1-2 GB download.

set -euo pipefail

ACC="GSE152495"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="${PROJECT_ROOT}/data_raw"
DATA_DIR="${PROJECT_ROOT}/data"

mkdir -p "$RAW_DIR" "$DATA_DIR"
cd "$RAW_DIR"

TAR="${ACC}_RAW.tar"
if [ ! -f "$TAR" ]; then
  echo ">>> Downloading ${ACC} supplementary archive..."
  wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/${ACC}/suppl/${ACC}_RAW.tar"
else
  echo ">>> ${TAR} already present, skipping download."
fi

echo ">>> Extracting..."
tar -xvf "$TAR"

echo
echo "=============================================================="
echo "Files extracted into ${RAW_DIR}:"
ls -1 "$RAW_DIR" | head -40
echo "=============================================================="
echo
echo "NEXT STEP — you must do this part by hand, once:"
echo
echo "1. Look at the filenames above. They will look something like"
echo "     GSM4616269_F1C_barcodes.tsv.gz"
echo "   The GSM number and the short code identify sex/treatment/replicate."
echo
echo "2. Open the GEO series page to see which GSM is which sample:"
echo "     https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${ACC}"
echo "   Cross-check against Supplemental Table S2 of paper.pdf."
echo
echo "3. Fill in the GSM -> sample mapping in tools/sample_map.tsv,"
echo "   then run:  bash tools/organize_data.sh"
echo
echo "DO NOT guess the mapping. Swapping male and female, or cocaine and"
echo "sucrose, produces a pipeline that runs perfectly and answers the"
echo "wrong question. This is the one step where a silent error costs you"
echo "the whole project."
