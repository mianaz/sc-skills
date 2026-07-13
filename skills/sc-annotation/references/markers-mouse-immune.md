# Canonical mouse immune/stromal markers (extensible)

Mouse symbols for the verification dotplot when markers.md (human) doesn't apply. Schema:
`cell_type | canonical_markers | citation`. Extend per project; keep a citation per row.
Work coarse compartments first, then fine subtypes (marker-verification.md).

| cell_type | canonical_markers | citation |
|---|---|---|
| Pan-immune | Ptprc (Cd45) | PMID 31178118 |
| T cell (pan) | Cd3d, Cd3e, Cd3g | PMID 31178118 |
| T CD4 | Cd4, Cd40lg, Cxcr3, Icos | PMID 31178118 |
| T CD8 | Cd8a, Cd8b1, Gzmk, Gzmb, Ifng | PMID 31178118 |
| Treg | Foxp3, Il2ra, Ctla4, Il7r | PMID 31178118 |
| T gamma-delta | Trdc, Trgv2, Rorc, Il23r, Il17a | PMID 28475889 |
| NK | Ncr1, Klrb1c, Klrk1, Il2rb, Nkg7 | PMID 31178118 |
| NKT | Nkg7, Klrb1c, Il4 | PMID 28475889 |
| ILC (pan) | Il7r, Id2, Gata3 (Lin-negative) | PMID 28475889 |
| ILC1 | Tbx21, Ncr1, Il7r, Il2rb, Cd160 | PMID 28475889 |
| ILC2 | Gata3, Areg, Il1rl1, Il5, Klrg1, Rora | PMID 28475889 |
| ILC3 | Rorc, Il7r, Csf2ra, Ltbr, Ncr1 | PMID 28475889 |
| B cell | Cd19, Ms4a1 (Cd20), Cd79a, Ighm, Ighd | PMID 31178118 |
| Plasma cell | Jchain, Mzb1, Xbp1, Ighg1 | PMID 31178118 |
| Monocyte | Ly6c2, Cd14, Ccr2, Fcgr1 (Cd64), Itgam | PMID 31178118 |
| Macrophage | Cd68, Adgre1 (F4/80), Mertk, Cd163, Mrc1, Csf1r | PMID 31178118 |
| Neutrophil | Ly6g, S100a8, Csf3r, Mpo, Cd177 | PMID 31178118 |
| Eosinophil | Siglecf, Ccr3, Il5ra, Prg2, Epx | PMID 28475889 |
| Basophil | Mcpt8, Fcer1a, Cd200r3, Il4 | PMID 28475889 |
| cDC1 | Xcr1, Batf3, Clec9a | PMID 27760337 |
| cDC2 | Itgax, Sirpa, Clec10a, Irf4 | PMID 27760337 |
| pDC | Siglech, Bst2, Tcf4, Irf7, Lilrb4a | PMID 27760337 |
| Fibroblast | Dcn, Col1a1, Col3a1, Sfrp4, Lum | PMID 36300619 |
| Endothelial | Pecam1, Vwf, Cldn5 | doi:10.1093/database/baz046 |
| Epithelial | Epcam, Krt8, Krt18 | PMID 36300619 |

Receptor↔symbol notes: Cd25=Il2ra, Cd127=Il7r, Cd11b=Itgam, Cd11c=Itgax, Cd122=Il2rb,
Nk1.1=Klrb1c, Nkg2d=Klrk1, Cd64=Fcgr1, F4/80=Adgre1, Cd20=Ms4a1.
