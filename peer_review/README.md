# Peer review — Milestone 2 and 3

## What goes in this folder

| File | When | Contents |
|---|---|---|
| `review_given.md` | Milestone 2 (by Sep 4) | The review you post on the assigned group's GitHub issue |
| `review_received.md` | Milestone 2 | Paste the review your group receives, so the rebuttal has a fixed reference |
| `rebuttal.md` | Milestone 3 (by Sep 9) | Point-by-point response, referencing commits |

---

## Review checklist

A review that only says "the figures look good" is worth nothing to the authors
and marks you down. Check things that could be *wrong*, and be specific about
where.

### Reproducibility
- [ ] Does the repo have a lockfile (`requirements.lock.txt` or `uv.lock`), or
      only unpinned package names?
- [ ] Does every script run from a clean clone, or does it depend on files that
      exist only on the author's machine?
- [ ] Are random seeds set? Leiden and UMAP are both stochastic — without a seed
      their "36 clusters" is not reproducible even by themselves.
- [ ] Is `data/` correctly excluded from git?

### Data handling
- [ ] Is the GSM → sample mapping documented and verified, not assumed? Ask how
      they confirmed sex and treatment labels.
- [ ] Are all 8 samples present, or did one silently fail to load?
- [ ] Does the cell count after QC bear a sensible relationship to the paper's
      86,224?

### Statistics — where most real errors hide
- [ ] **Is DE computed on log-normalized data, or accidentally on scaled
      z-scores?** Check whether `use_raw=True` is set, or whether `adata.raw`
      was assigned before `sc.pp.scale`. This is the single most common bug and
      it silently corrupts every fold change.
- [ ] Is `reference=` set correctly, so positive log2FC means up in cocaine?
      If reversed, every biological statement in the Discussion inverts.
- [ ] **Fold-change units:** scanpy reports log2, the paper reports natural log.
      Did they compare |log2FC| > 1 against the paper's |ln FC| > 1 as if the
      two were the same number? They are not — the paper's cut equals
      |log2FC| > 1.44.
- [ ] Is multiple-testing correction applied, and is the method stated?
      (Paper: Bonferroni. Scanpy default: Benjamini–Hochberg. Different results.)
- [ ] Is pseudoreplication acknowledged? Cell-level Wilcoxon treats ~80,000
      cells as independent replicates when there are really 2 per group.

### Biological interpretation
- [ ] **Are clusters matched to the paper by marker genes, or by number?**
      Their C11 is Kenyon cells; the reviewee's cluster 11 is almost certainly
      something else. Matching by number is a real error that looks like a
      result.
- [ ] Is every annotation supported by a named marker gene?
- [ ] Are unannotated clusters labelled honestly, or given a confident guess?
- [ ] Do the Enrichr libraries actually exist for fly? `KEGG_2021_Human` does
      not, and returns nothing useful without erroring.

### Reporting
- [ ] Are the differences from the paper (Seurat vs Scanpy, SCTransform vs
      CPM+log1p, MAST vs Wilcoxon) acknowledged, or is a mismatch presented as
      a failure?
- [ ] Is the male:female DE ratio reported, not just raw counts? The ratio is
      the reproducible claim.
- [ ] Is the AI Usage Disclosure appendix present and specific — actual prompts,
      not "we used ChatGPT for help"?

---

## Tone

Separate three things and say which is which:

1. **Errors** — something is wrong and the conclusion does not follow.
2. **Concerns** — the result may hold, but the evidence shown does not establish it.
3. **Suggestions** — it would be better with X, but it is not wrong without it.

Anchor each point to a file and line or a figure number. "The clustering seems
off" is not actionable; "`03_cluster.py` line 47 calls `sc.tl.leiden` without
`random_state`, so the 36 clusters in Figure 2 are not reproducible" is.

---

## Rebuttal format (Milestone 3)

For each point received:

> **Reviewer comment:** *(quote it verbatim)*
>
> **Response:** Accepted / Partially accepted / Respectfully disagree
>
> **Action taken:** what changed, and the commit hash.
>
> **If disagreeing:** the reason, with evidence. Disagreeing is legitimate when
> you are right — but the burden is on you to show it, not to restate your
> original claim.
