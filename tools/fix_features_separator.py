#!/usr/bin/env python3
"""
fix_features_separator.py — convert space-delimited 10x metadata files to TSV.

THE PROBLEM
`sc.read_10x_mtx` parses features.tsv.gz with sep="\\t", per the CellRanger
specification. The files supplied for this project are SPACE-delimited, so every
field lands in column 0 and scanpy fails with `KeyError: 1` when it reaches
`genes[1]` looking for the gene symbol.

Diagnosis that identified it:
    zcat data/Female_Cocaine_1/features.tsv.gz | head -1 | awk -F'\\t' '{print NF}'
    -> 1          (should be 3)

WHY CONVERT THE FILES RATHER THAN PATCH THE READER
Writing a bespoke loader means every downstream tool that expects standard 10x
input needs the same workaround, and a peer reviewer cloning the repo hits the
same wall. Normalizing the input once keeps the pipeline standard.

The third column is "Gene Expression" -- which itself contains a space. So we
split on whitespace with maxsplit=2, giving exactly [gene_id, symbol, type].
Splitting naively on every space would produce four fields and corrupt the
feature type.

SAFETY
Originals are backed up to <sample>/original/ before anything is written, and
the script refuses to run twice on the same file. Re-running is harmless.

Run:  python tools/fix_features_separator.py           # report only
      python tools/fix_features_separator.py --apply   # convert
"""

import argparse
import gzip
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from config import DATA_DIR


def read_first_line(path):
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n")


def needs_fixing(path):
    """True if the file has no tabs but does have whitespace-separated fields."""
    first = read_first_line(path)
    if "\t" in first:
        return False
    return len(first.split()) >= 2


def convert(path, backup_dir, dry_run=True):
    """Rewrite a space-delimited file as tab-delimited, preserving field count."""
    first = read_first_line(path)
    n_fields = len(first.split(None, 2))

    print(f"    fields detected: {n_fields}")
    print(f"    first line: {first[:70]}")

    if dry_run:
        fixed = "\t".join(first.split(None, 2))
        print(f"    would become: {fixed[:70]}")
        return

    backup_dir.mkdir(exist_ok=True)
    backup = backup_dir / path.name
    if not backup.exists():
        shutil.copy2(path, backup)
        print(f"    backed up -> {backup}")

    lines = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            # maxsplit=2 keeps "Gene Expression" intact as one field
            lines.append("\t".join(line.split(None, 2)))

    tmp = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(tmp, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    tmp.replace(path)
    print(f"    converted {len(lines):,} lines")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="actually convert (default is report only)")
    args = ap.parse_args()

    samples = sorted(p for p in DATA_DIR.iterdir() if p.is_dir() and p.name != "original")
    if not samples:
        sys.exit(f"No sample folders under {DATA_DIR}")

    to_fix = 0
    for sample in samples:
        print(f"\n{sample.name}")
        for fname in ["features.tsv.gz", "barcodes.tsv.gz"]:
            path = sample / fname
            if not path.exists():
                print(f"  {fname}: MISSING")
                continue

            if needs_fixing(path):
                print(f"  {fname}: space-delimited, needs conversion")
                convert(path, sample / "original", dry_run=not args.apply)
                to_fix += 1
            else:
                first = read_first_line(path)
                n = len(first.split("\t"))
                # barcodes.tsv is single-column by design -- nothing to fix
                print(f"  {fname}: OK ({n} tab-separated field"
                      f"{'s' if n != 1 else ''})")

    print()
    if not args.apply and to_fix:
        print(f"{to_fix} file(s) need conversion. Re-run with --apply to do it.")
        print("Originals are backed up to data/<sample>/original/ first.")
    elif args.apply:
        print(f"Converted {to_fix} file(s).")
        print("\nIMPORTANT: delete any stale scanpy cache before re-reading:")
        print("    rm -rf cache/")
        print("\nRecord this conversion in docs/data_provenance.md -- you have")
        print("modified the supplied data, and that belongs in the record.")
    else:
        print("Nothing to fix; all files are already tab-separated.")


if __name__ == "__main__":
    main()
