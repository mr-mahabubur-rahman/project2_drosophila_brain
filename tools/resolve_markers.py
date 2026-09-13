#!/usr/bin/env python3
"""
resolve_markers.py — check every marker in config.CANONICAL_MARKERS against the
symbols actually used in this dataset's reference annotation.

THE PROBLEM
The CellRanger reference used for this dataset does not use current FlyBase
symbols throughout. For a subset of genes it substitutes vertebrate ortholog
names:

    FlyBase ID   fly symbol   this reference calls it
    CG9887       VGlut        VGlut1
    CG9533       rut          Adcy1
    CG32498      dnc          Pde4
    CG9122       Trh          Trhn
    CG30035      Tret1-1      Tret1

A marker that is silently absent gives you a cluster you cannot label. Far
worse, a marker that silently resolves to the WRONG GENE gives you a label you
will defend in your report.

THE TRAP THIS DATASET CONTAINS
`trh` (lowercase) in this reference is CG42865 = trachealess, a tracheal
transcription factor. `Trh` (capital) is CG9122 = tryptophan hydroxylase, the
serotonin synthesis enzyme. A case-insensitive match finds trachealess and
reports success. Drosophila nomenclature is case-sensitive and this is not an
edge case -- it is the single most likely way to mislabel a cluster here.

WHAT THIS SCRIPT DOES
For each configured marker it reports: exact symbol match, case-insensitive
match (flagged as suspect), alias resolution by FlyBase ID, or not found. It
writes a worksheet for you to verify each resolution on FlyBase.

IT DOES NOT AUTO-EDIT YOUR CONFIG. Every substitution must be checked by you
against flybase.org before you rely on it. That is the whole lesson of the
trh case.

Run:  python tools/resolve_markers.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import pandas as pd

from config import DATA_DIR, TABLE_DIR, CANONICAL_MARKERS

# CONFIRMED by direct lookup in this dataset's features.tsv.
# fly symbol -> (FlyBase CG id, symbol used by this reference)
CONFIRMED_ALIASES = {
    "VGlut":   ("CG9887",  "VGlut1"),
    "rut":     ("CG9533",  "Adcy1"),
    "dnc":     ("CG32498", "Pde4"),
    "Trh":     ("CG9122",  "Trhn"),
    "Tret1-1": ("CG30035", "Tret1"),
}

# Symbols that resolve to a DIFFERENT gene in this reference. Never use these.
DANGEROUS = {
    "trh": "CG42865 = trachealess, a tracheal TF. NOT tryptophan hydroxylase "
           "(that is Trh / CG9122, called 'Trhn' here).",
}


def load_features():
    sample = next(p for p in sorted(DATA_DIR.iterdir())
                  if p.is_dir() and (p / "features.tsv.gz").exists())
    feat = pd.read_csv(
        sample / "features.tsv.gz", sep="\t", header=None,
        names=["gene_id", "symbol", "type"], keep_default_na=False,
    )
    print(f"Reference annotation from {sample.name}: {len(feat):,} features\n")
    return feat


def resolve(feat, marker):
    """Return (status, resolved_symbol, gene_id, note)."""
    exact = feat[feat["symbol"] == marker]
    if len(exact):
        r = exact.iloc[0]
        return "OK", marker, r["gene_id"], ""

    if marker in CONFIRMED_ALIASES:
        cg, alias = CONFIRMED_ALIASES[marker]
        hit = feat[feat["symbol"] == alias]
        if len(hit):
            return ("ALIAS", alias, hit.iloc[0]["gene_id"],
                    f"{marker} ({cg}) is named {alias} in this reference")
        return "MISSING", None, None, f"alias {alias} also absent"

    ci = feat[feat["symbol"].str.lower() == marker.lower()]
    if len(ci):
        r = ci.iloc[0]
        note = DANGEROUS.get(r["symbol"], "")
        if note:
            return "DANGER", r["symbol"], r["gene_id"], note
        return ("CASE", r["symbol"], r["gene_id"],
                f"only a case-insensitive match -- VERIFY on FlyBase that "
                f"{r['symbol']} is the same gene as {marker}")

    similar = feat[feat["symbol"].str.contains(marker, case=False, regex=False)]
    cands = similar["symbol"].tolist()[:5]
    return "MISSING", None, None, f"candidates: {cands}" if cands else "no similar symbols"


def main():
    feat = load_features()
    rows = []

    for cell_type, markers in CANONICAL_MARKERS.items():
        print(f"{cell_type}")
        for m in markers:
            status, sym, gid, note = resolve(feat, m)
            flag = {"OK": "  ok   ", "ALIAS": " ALIAS ",
                    "CASE": " CASE? ", "DANGER": " DANGER",
                    "MISSING": "MISSING"}[status]
            line = f"  {flag} {m:<12}"
            if sym and sym != m:
                line += f" -> {sym:<12} ({gid})"
            elif sym:
                line += f"    {gid}"
            print(line)
            if note:
                print(f"          {note}")
            rows.append({
                "cell_type": cell_type, "configured_marker": m,
                "status": status, "symbol_in_reference": sym,
                "gene_id": gid, "note": note,
                "VERIFIED_ON_FLYBASE": "",
            })
        print()

    df = pd.DataFrame(rows)
    out = TABLE_DIR / "marker_resolution.csv"
    df.to_csv(out, index=False)

    # ---- Summary ---------------------------------------------------------
    counts = df["status"].value_counts()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for status in ["OK", "ALIAS", "CASE", "DANGER", "MISSING"]:
        if status in counts:
            print(f"  {status:<8} {counts[status]:>3}")

    danger = df[df["status"] == "DANGER"]
    if len(danger):
        print("\nDANGER -- these resolve to a DIFFERENT gene. Remove them:")
        for _, r in danger.iterrows():
            print(f"  {r['configured_marker']}: {r['note']}")

    case = df[df["status"] == "CASE"]
    if len(case):
        print("\nCASE-ONLY matches -- verify each on FlyBase before trusting:")
        for _, r in case.iterrows():
            print(f"  {r['configured_marker']} -> {r['symbol_in_reference']} "
                  f"({r['gene_id']})")

    missing = df[df["status"] == "MISSING"]
    if len(missing):
        print("\nNOT FOUND -- look each up on FlyBase, find its CG number, then")
        print("search this reference for that ID to get the symbol it uses:")
        for _, r in missing.iterrows():
            print(f"  {r['configured_marker']:<12} {r['note']}")

    print(f"\nWorksheet: {out}")
    print("\nFill in VERIFIED_ON_FLYBASE for every ALIAS, CASE and MISSING row,")
    print("then update config.CANONICAL_MARKERS. Do not skip the verification --")
    print("the trh case shows what a plausible-looking wrong match costs.")


if __name__ == "__main__":
    main()
