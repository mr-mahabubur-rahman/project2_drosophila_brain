#!/usr/bin/env bash
# organize_data.sh — turn flat GEO files into per-sample 10x folders.
#
# Reads tools/sample_map.tsv (two columns: GSM_prefix <TAB> Sample_Name)
# and builds data/<Sample_Name>/{barcodes,features,matrix} for each row.
#
# Usage:  bash tools/organize_data.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="${PROJECT_ROOT}/data_raw"
DATA_DIR="${PROJECT_ROOT}/data"
MAP="${PROJECT_ROOT}/tools/sample_map.tsv"

[ -f "$MAP" ] || { echo "ERROR: $MAP not found."; exit 1; }

while IFS=$'\t' read -r gsm sample; do
  # skip comments and blanks
  [[ -z "${gsm// }" || "$gsm" == \#* ]] && continue

  echo ">>> ${gsm} -> ${sample}"
  mkdir -p "${DATA_DIR}/${sample}"

  for kind in barcodes features matrix; do
    # Match e.g. GSM4616269_F1C_barcodes.tsv.gz  (also handles 'genes' naming)
    src=$(find "$RAW_DIR" -maxdepth 1 -name "${gsm}*${kind}*" | head -1)

    if [ -z "$src" ] && [ "$kind" = "features" ]; then
      # Older CellRanger (v2) wrote genes.tsv.gz instead of features.tsv.gz
      src=$(find "$RAW_DIR" -maxdepth 1 -name "${gsm}*genes*" | head -1)
    fi

    if [ -z "$src" ]; then
      echo "    WARNING: no '${kind}' file found for ${gsm}"
      continue
    fi

    ext="tsv.gz"
    [ "$kind" = "matrix" ] && ext="mtx.gz"
    cp "$src" "${DATA_DIR}/${sample}/${kind}.${ext}"
    echo "    $(basename "$src") -> ${kind}.${ext}"
  done
done < "$MAP"

echo
echo ">>> Verifying..."
ok=0
for d in "${DATA_DIR}"/*/; do
  name=$(basename "$d")
  if [ -f "${d}barcodes.tsv.gz" ] && [ -f "${d}features.tsv.gz" ] && [ -f "${d}matrix.mtx.gz" ]; then
    echo "  OK      ${name}"
    ok=$((ok+1))
  else
    echo "  MISSING ${name} — incomplete, step 01 will fail"
  fi
done

echo
echo "${ok} complete samples (expected 8)."
[ "$ok" -eq 8 ] && echo "Ready. Run: uv run scripts/01_load_data.py"
