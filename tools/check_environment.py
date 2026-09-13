#!/usr/bin/env python3
"""
check_environment.py — verify the environment and project layout before you
burn an hour on a run that was going to fail at minute 55.

Checks, in order: Python version, every required package, the Leiden backend
(the most common silent failure), the folder layout, the reference PDFs, the
data folders, and available RAM against what step 03 will need.

Run:  uv run tools/check_environment.py
"""

import importlib
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

OK, WARN, FAIL = "  OK   ", "  WARN ", "  FAIL "
problems, warnings_ = [], []


def check(label, condition, detail="", fatal=True):
    if condition:
        print(f"{OK} {label} {detail}")
    else:
        print(f"{FAIL if fatal else WARN} {label} {detail}")
        (problems if fatal else warnings_).append(label)
    return condition


print("=" * 70)
print("PYTHON")
print("=" * 70)
v = sys.version_info
check(f"Python {v.major}.{v.minor}.{v.micro}", v >= (3, 9),
      "(project specifies 3.11)")
if v[:2] != (3, 11):
    warnings_.append("Python is not 3.11")
    print(f"{WARN} not 3.11 — peer reviewers reproducing your lockfile may differ")

print()
print("=" * 70)
print("PACKAGES")
print("=" * 70)
required = ["scanpy", "anndata", "numpy", "pandas", "scipy", "matplotlib", "gseapy"]
optional = ["seaborn", "harmonypy", "jupyter"]

for pkg in required:
    try:
        m = importlib.import_module(pkg)
        check(pkg, True, f"v{getattr(m, '__version__', '?')}")
    except ImportError:
        check(pkg, False, "— uv pip install -r requirements.txt")

for pkg in optional:
    try:
        m = importlib.import_module(pkg)
        print(f"{OK} {pkg} v{getattr(m, '__version__', '?')} (optional)")
    except ImportError:
        print(f"{WARN} {pkg} missing (optional)")

print()
print("=" * 70)
print("LEIDEN BACKEND")
print("=" * 70)
print("scanpy does NOT install these itself. sc.tl.leiden fails at step 03")
print("with an unhelpful error if they are absent.\n")
for pkg in ["igraph", "leidenalg"]:
    try:
        m = importlib.import_module(pkg)
        check(pkg, True, f"v{getattr(m, '__version__', '?')}")
    except ImportError:
        check(pkg, False, f"— uv pip install {pkg}")

print()
print("=" * 70)
print("PROJECT LAYOUT")
print("=" * 70)
try:
    import config
    root = config.PROJECT_ROOT
    print(f"  Root: {root}\n")
    for name in ["data", "docs", "notebooks", "peer_review",
                 "reports", "results", "scripts", "tools"]:
        check(f"{name}/", (root / name).is_dir(), fatal=(name != "data"))

    print()
    check("docs/paper.pdf", config.PAPER_PDF.exists(), fatal=False)
    check("docs/project2_guide.pdf", config.GUIDE_PDF.exists(), fatal=False)
except Exception as e:
    check("config.py importable", False, f"— {e}")
    config = None

print()
print("=" * 70)
print("DATA")
print("=" * 70)
if config is not None and config.DATA_DIR.exists():
    samples = sorted(p for p in config.DATA_DIR.iterdir() if p.is_dir())
    expected = {f"{s}_{t}_{r}" for s in ("Female", "Male")
                for t in ("Cocaine", "Sucrose") for r in ("1", "2")}
    found = set()
    for p in samples:
        complete = all((p / f).exists() for f in
                       ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"])
        print(f"{OK if complete else WARN} {p.name}"
              f"{'' if complete else ' — missing one of the 10x triplet'}")
        if complete:
            found.add(p.name)

    missing = expected - found
    if missing:
        print(f"\n{WARN} Missing samples: {sorted(missing)}")
        print("        Run: bash tools/download_data.sh, fill tools/sample_map.tsv,")
        print("        then bash tools/organize_data.sh")
    else:
        print(f"\n{OK} All 8 samples present and complete.")
else:
    print(f"{WARN} data/ not found or empty — this is expected before download.")

print()
print("=" * 70)
print("RESOURCES")
print("=" * 70)
try:
    import os
    total_ram = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9
    is_wsl = "microsoft" in Path("/proc/version").read_text().lower() \
        if Path("/proc/version").exists() else False

    print(f"  RAM visible to Linux: {total_ram:.1f} GB{'  (WSL2)' if is_wsl else ''}")

    if is_wsl:
        print("\n  NOTE: on WSL2 this is what the Linux VM was GIVEN, not your host RAM.")
        print("  WSL2 defaults to half the host (or 8 GB, whichever is less).")
        print("  On a 16 GB PC that means Linux sees 8 GB unless you override it.")
        print("  Fix: create C:\\Users\\<you>\\.wslconfig containing")
        print("       [wsl2]")
        print("       memory=12GB")
        print("       swap=8GB")
        print("  then run  wsl --shutdown  in PowerShell and reopen your terminal.\n")

    if total_ram < 8:
        print(f"{FAIL} Under 8 GB. Step 03 will be killed by the OOM killer.")
        print("        Raise the WSL allocation (above), then set REGRESS_OUT=False.")
        problems.append("insufficient RAM")
    elif total_ram < 11:
        print(f"{WARN} 8–11 GB. Step 03 may survive with REGRESS_OUT=False, but")
        print("        raising the WSL allocation to 12 GB is the better fix.")
        warnings_.append("tight RAM")
    else:
        print(f"{OK} Enough for the full ~80k-cell run with REGRESS_OUT=True.")
        print("        Close Chrome and other heavy applications before step 03 anyway.")
except Exception as e:
    print(f"{WARN} Could not read RAM ({e}).")

free_gb = shutil.disk_usage(Path.home()).free / 1e9
print(f"  Free disk: {free_gb:.1f} GB")
if free_gb < 20:
    print(f"{WARN} Under 20 GB. Raw data (~2 GB) plus four .h5ad checkpoints")
    print("        can approach 15 GB.")

if str(Path.cwd()).startswith("/mnt/c") or str(Path.cwd()).startswith("/mnt/d"):
    print(f"\n{WARN} You are on the Windows filesystem via WSL.")
    print("        Reading thousands of matrix entries across the /mnt/ boundary is")
    print("        roughly 10x slower. Move the project into ~/ instead.")
    warnings_.append("running from /mnt/c")

print()
print("=" * 70)
if problems:
    print(f"BLOCKED — fix these first: {', '.join(problems)}")
    sys.exit(1)
elif warnings_:
    print(f"READY, with caveats: {', '.join(warnings_)}")
else:
    print("READY. Next: uv run scripts/01_load_data.py")
print("=" * 70)
