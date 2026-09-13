# Cluster annotation — reference and draft

## The paper's 36 clusters

From `new.cluster.ids` in the authors' `Rcode_for_analyses.R`, cross-referenced
with the cell-type labels in Figure 2 of `docs/paper.pdf`. The triplet is the
top-3 marker genes they used to define each cluster.

| Paper | Marker triplet | Cell type (Fig. 2) |
|---|---|---|
| C0 | VGlut/CG2269/CG32017 | Glutamatergic neurons type 1 |
| C1 | Gad1/CG14989/CG32017 | Central brain GABAergic neurons type 1 |
| C2 | spab/jeb/CG31221 | Neuropeptide/cholinergic neurons |
| C3 | VGlut/CG34355/CG9650 | Glutamatergic |
| C4 | SoxN/CG9650/klg | Optic lobe type 1 |
| C5 | pros/dati/spab | Central brain A type 1 |
| C6 | toy/bi | — |
| C7 | pros/br/beat-IIIC | Central brain A type 2 |
| C8 | sosie/acj6/nAChRalpha6 | Olfactory projection neurons |
| C9 | Gad1/Lim3/CG14989 | GABAergic neurons type 2 |
| C10 | mbl/Imp/CG31345 | Central brain B cholinergic neurons |
| **C11** | **jdp/Pka-R2/Rgk1** | **Kenyon cells type 1** |
| C12 | ct/Gad1/CG14989 | GABAergic neurons type 3 |
| C13 | Obp44a/CG8369/CG15201 | Glia type 1 |
| C14 | pdm3/vvl/CG18598 | Optic lobe type 2 |
| C15 | Pka-C3/sosie/CG2016 | — |
| **C16** | **CNMaR/SoxN/CG42750** | **Unannotated (antennal/optic lobe mix)** |
| **C17** | **Eaat1/CG2016/CG15522** | **Astrocytes** |
| C18 | DIP-theta/aop/CG42458 | Optic lobe type 3 |
| **C19** | **Gs2/CG1552/CG8369** | **Glial type 2** |
| C20 | Pka-R2/crb/Rgk1 | Kenyon cells type 2 |
| C21 | VGlut/mbl/Proc | — |
| **C22** | **Cys/SPARC/CG3168** | **Surface glia and fat body** |
| C23 | Nos/Octbeta1R/CG14989 | GABAergic and monoaminergic |
| C24 | Vmat/Hsp23/CG4577 | Central brain type 2 monoaminergic |
| C25 | bi/shakB/Octbeta1R | Glutamatergic neurons type 2 |
| C26 | shakB/CG3134/CG18598 | — |
| C27 | Dh31/Nep1/CG14687 | Mushroom body and ellipsoid body |
| C28 | Frq1/CG43795/CG14274 | — |
| C29 | DAT/Vmat/ple | Dopaminergic neurons |
| C30 | Arr2/Arr1/ninaE | Photoreceptor cells |
| C31 | CG2016/CG6044/CG42540 | — |
| C32 | kn/sNPF/lncR:CR45223 | — |
| C33 | VGlut/Ca-alpha1T/CG8034 | Glutamatergic neurons type 3 |
| C34 | Tk/Nplp1/nAChRalpha6 | Tachykinin/neuropeptidergic |
| C35 | sNPF/CG14989/lncR:CR45223 | Central brain GABAergic |

**Bold rows are the paper's four strongest cocaine responders.** Finding your
equivalents of C11, C16, C17 and C22 is what Research Question 2 turns on.

Note: the paper's C11 (Kenyon cells) is defined by **jdp/Pka-R2/Rgk1**, not by
the `ey`/`Fas2`/`rut`/`dnc` panel in the project guide. Both are legitimate
Kenyon-cell markers, but matching *their* cluster requires *their* genes.

---

## Draft assignments for your 30 clusters

**These are suggestions, not answers.** Each needs checking against the dotplots
in `results/figures/` before you write it into the worksheet. Confidence is my
assessment of how much checking each needs.

| You | n | Your top-3 | Draft assignment | Paper match | Confidence |
|---|---|---|---|---|---|
| 0 | 8,502 | lncRNA:mimi, toy, VAChT | Cholinergic neurons | C6 (toy/bi) | medium |
| 1 | 3,196 | CG8369, lncRNA:CR34335, MtnA | Glia type 1 | C13 (Obp44a/CG8369/…) | medium |
| 2 | 143 | CG9394, CG1552, Gs2 | Glial type 2 | C19 (Gs2/CG1552/CG8369) | **high** |
| 3 | 1,047 | vvl, nAChRalpha6, CG42541 | Optic lobe type 2 | C14 (pdm3/vvl/CG18598) | **high** |
| 4 | 7,543 | VGlut1, cpx, CG34355 | Glutamatergic type 3 | C3 (VGlut/CG34355/CG9650) | **high** |
| 5 | 695 | lncRNA:CR34335, CG13631, eEF1alpha1 | ? high-expression cluster | — | low |
| 6 | 10,273 | VGlut1, IA-2, CG2269 | Glutamatergic type 1 | C0 (VGlut/CG2269/CG32017) | **high** |
| 7 | 9,058 | CG31221, IA-2, lncRNA:Hsromega | Neuropeptide/cholinergic | C2 (spab/jeb/CG31221) | medium |
| 8 | 1,913 | Vmat, ple, DAT | **Dopaminergic neurons** | C29 (DAT/Vmat/ple) | **very high** |
| 9 | 1,775 | disco, disco-r, CG31221 | Optic lobe | — | medium |
| 10 | 1,462 | nAChRalpha6, sosie, beat-Ic | Olfactory projection neurons | C8 (sosie/acj6/nAChRα6) | **high** |
| 11 | 3,032 | Imp, futsch, sesB | Central brain B cholinergic | C10 (mbl/Imp/CG31345) | **high** |
| 12 | 3,822 | mub, jdp, Pka-C1 | **Kenyon cells type 1** | C11 (jdp/Pka-R2/Rgk1) | medium — **verify** |
| 13 | 9,472 | CG14989, Gad1, CG10804 | Central brain GABAergic type 1 | C1 (Gad1/CG14989/…) | **high** |
| 14 | 893 | dths, cpx, CG18598 | ? | C26 (shakB/CG3134/CG18598) | low |
| 15 | 2,221 | SoxN, CG42750, acj6 | **Unannotated (antennal/optic)** | C16 (CNMaR/SoxN/CG42750) | **high** |
| 16 | 7,093 | pros, VAChT, acj6 | Central brain A type 1 | C5 (pros/dati/spab) | medium |
| 17 | 3,408 | ct, Rbp6, Jhbp2 | GABAergic type 3 | C12 (ct/Gad1/CG14989) | medium |
| 18 | 1,998 | Rbp6, Cbp53E, ct | GABAergic type 3 (2nd) | C12 | low |
| 19 | 526 | CG1552, nrv2, CG8369 | Glial type 2 (2nd) | C19 | medium |
| 20 | 1,542 | CG14989, Octbeta1R, Bx | GABAergic and monoaminergic | C23 (Nos/Octbeta1R/…) | medium |
| 21 | 474 | Pdh, fabp, pkm | Pigment/retinal support | — | low |
| 22 | 1,181 | beat-IIa, Con, nAChRalpha7 | ? | — | low |
| 23 | 526 | CG43795, CG42342, hth | ? | C28 (Frq1/CG43795/CG14274) | low |
| 24 | 1,915 | Jhbp2, Eaat1, CG15522 | **Astrocytes** | C17 (Eaat1/CG2016/CG15522) | **very high** |
| 25 | 822 | lncRNA:CR34335, Oda, CG3168 | **Surface glia / fat body** | C22 (Cys/SPARC/CG3168) | **high** |
| 26 | 444 | Arr2, Arr1, ninaE | **Photoreceptor cells** | C30 (Arr2/Arr1/ninaE) | **very high** |
| 27 | 464 | VGlut1, Argk1, Fim | Glutamatergic (small) | — | low |
| 28 | 427 | jdp, 14-3-3zeta, Tomosyn | Kenyon cells (2nd) | C11/C20 | low — **verify** |
| 29 | 310 | VGlut1, CG43689, Lmx1a | Glutamatergic (small) | — | low |

---

## The four that matter most for RQ2

Your likely equivalents of the paper's top responders:

| Paper | Cell type | Your cluster | Basis |
|---|---|---|---|
| C11 | Kenyon cells | **12** (maybe 28) | `jdp` shared; check `Pka-R2`, `Rgk1` |
| C16 | Unannotated | **15** | `SoxN` + `CG42750` both shared |
| C17 | Astrocytes | **24** | `Eaat1` + `CG15522` both shared |
| C22 | Surface glia | **25** | `CG3168` shared; check `Cys`, `SPARC` |

Three of these are solid on two shared markers each. **Cluster 12 as Kenyon
cells rests on one gene (`jdp`) and needs confirming.** Check `Pka-R2`, `Rgk1`,
`crb` in clusters 12 and 28 before committing — this is the single assignment
your RQ2 conclusions depend on most.

---

## Verification to run before filling the worksheet

```python
import scanpy as sc, pandas as pd, sys
sys.path.insert(0, 'scripts'); import config

adata = sc.read_h5ad(config.H5AD_ANNOTATED)
key = config.LEIDEN_KEY

# Paper's defining genes for the four key clusters
panels = {
    'C11 Kenyon':      ['jdp', 'Pka-R2', 'Rgk1', 'crb'],
    'C16 unannotated': ['CNMaR', 'SoxN', 'CG42750'],
    'C17 astrocytes':  ['Eaat1', 'CG2016', 'CG15522'],
    'C22 surface glia':['Cys', 'SPARC', 'CG3168'],
}
genes = [g for gs in panels.values() for g in gs if g in adata.raw.var_names]
missing = [g for gs in panels.values() for g in gs if g not in adata.raw.var_names]
print('not in reference:', missing)

df = sc.get.obs_df(adata, keys=genes, use_raw=True)
df[key] = adata.obs[key].values
m = df.groupby(key, observed=True).mean().round(3)

for label, panel in panels.items():
    cols = [g for g in panel if g in m.columns]
    if not cols:
        continue
    print(f'\n{label} — top 5 clusters by mean of {cols}:')
    print(m[cols].mean(axis=1).sort_values(ascending=False).head(5).round(3).to_string())
```

If cluster 12 tops the C11 panel, the Kenyon-cell assignment is confirmed. If
cluster 28 or something else tops it, reassign.

---

## Filling the worksheet

1. Open `results/tables/annotation_worksheet.csv`
2. For each cluster, fill `cell_type_ANNOTATE_ME`, `paper_cluster_match`, and
   `evidence_notes` — name the genes you based it on
3. Label clusters you cannot resolve as `Unannotated`. The paper did this for
   their C16 and said so; an honest blank beats a confident guess
4. Save as `annotation_worksheet_FILLED.csv` in the same folder
5. Re-run `python scripts/04_annotate_clusters.py`

Remember the symbol aliases when reading the dotplots: `VGlut1` = VGlut,
`Adcy1` = rut, `Pde4` = dnc, `Trhn` = Trh, `Tret1` = Tret1-1.
