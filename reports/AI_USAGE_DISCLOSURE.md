# Appendix: AI Usage Disclosure

*Required by the Project 2 academic integrity policy. Fill this in as you work.*

## Tools used

| Tool | Model | Used for |
|---|---|---|
| Claude (Anthropic) | Claude Opus 5 | Pipeline scaffolding, parameter rationale, debugging |
|  |  |  |

## Prompt log

| # | Date | Prompt (verbatim) | Purpose | What I changed / verified |
|---|---|---|---|---|
| 1 | 2026-09-12 | "here is the group to paper and guideline. Now prepare the whole pipeline from environment setup. Do you need anything else?" | Generate the initial six-script Scanpy pipeline and environment setup from the supplied paper and project guide | *(record here: which parameters you kept, which you changed, and what you checked)* |
| 2 |  |  |  |  |

## Validation notes

The policy states you are solely responsible for validating all outputs. Record
your checks here. Suggested minimum:

- [ ] **Sample mapping verified** against the GEO series page and Supplemental
      Table S2, plus the `roX1`/`roX2` sex check (see README §2).
- [ ] **Marker gene names verified** against FlyBase. AI-suggested marker lists
      are a known hallucination risk — `alrm`, `wrapper`, `Tret1-1`, `Mdr65`,
      `zyd` and `Gat` in `config.py` were added by the AI beyond the guide's
      table and must each be confirmed.
- [ ] **Cluster annotations are my own**, derived from the marker dotplots and
      the paper's Supplemental Table S4 — not from an AI-suggested label.
- [ ] **Enrichr library names checked** against `gp.get_library_name(organism='fly')`
      (the project guide's `KEGG_2021_Human` does not exist in the Fly catalogue).
- [ ] **Fold-change units checked**: scanpy reports log2; the paper reports
      natural log. Confirmed the conversion used in `config.py`.
- [ ] **DE direction checked**: `reference='Sucrose'` means positive log2FC =
      up in cocaine. Spot-verified on one gene by hand.
- [ ] **Every parameter in `config.py` reviewed** and I can justify each at
      grading.
- [ ] **Statistical caveat understood**: cell-level Wilcoxon treats cells as
      independent replicates (pseudoreplication); `05b` cross-checks against
      pseudobulk.

## Statement

I understand all code submitted in this project and can explain every line and
parameter choice. Where AI assistance was used, I have verified the outputs
against primary sources (FlyBase, the GEO record, and Baker et al. 2021).

Signed: ______________________  Date: ____________
