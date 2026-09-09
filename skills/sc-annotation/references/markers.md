# Canonical marker reference (extensible)

Schema: `cell_type | canonical_markers | tissue/context | citation (PMID/DOI)`.
Seeded from CellMarker 2.0 / PanglaoDB / Azimuth references. Extend per project; keep a citation for every row.

| cell_type | canonical_markers | tissue/context | citation |
|---|---|---|---|
| T cell (CD4) | CD3D, CD3E, CD4, IL7R | PBMC/general | PMID 34062119 |
| T cell (CD8) | CD3D, CD8A, CD8B, GZMK | PBMC/general | PMID 34062119 |
| NK cell | NCAM1, NKG7, GNLY, KLRD1 | PBMC/general | PMID 34062119 |
| B cell | CD19, MS4A1, CD79A, CD79B | PBMC/general | PMID 34062119 |
| Plasma cell | SDC1, MZB1, XBP1, PRDM1 | general | PMID 36300619 |
| Monocyte/Macrophage | CD14, LYZ, FCGR3A, CD68 | PBMC/tissue | PMID 34062119 |
| Dendritic cell | FCER1A, CST3, CLEC9A, LILRA4 | PBMC/tissue | PMID 34062119 |
| Mast cell | TPSAB1, CPA3, KIT | tissue | doi:10.1093/database/baz046 |
| Endothelial | PECAM1, VWF, CLDN5 | tissue | doi:10.1093/database/baz046 |
| Fibroblast | COL1A1, DCN, LUM, PDGFRB | tissue | doi:10.1093/database/baz046 |
| Epithelial | EPCAM, KRT8, KRT18 | tissue | PMID 36300619 |

## Liver / tumour-immune-microenvironment (HCC, ICC)

Vetted lineage panels from a >1M-cell liver-cancer atlas. Use alongside the general rows above.
Citation for this block: Xue et al. *Nature* 2022, doi:10.1038/s41586-022-05400-x
(source: sc-paper-distill/papers/2022-scPLC.md F2).

| cell_type | canonical_markers | tissue/context | citation |
|---|---|---|---|
| Neutrophil / TAN | CSF3R, S100A8, S100A9, FCGR3B, CXCR2 | liver tumour | doi:10.1038/s41586-022-05400-x |
| Macrophage (incl. Kupffer) | CD68, CD163, C1QC, MARCO; Kupffer-resident VSIG4, CD5L, TIMD4 | liver | doi:10.1038/s41586-022-05400-x |
| Monocyte | FCN1, S100A8, CD14; CD163 (mono-mac) | liver | doi:10.1038/s41586-022-05400-x |
| Dendritic cell | LILRA4 (pDC), CLEC9A (cDC1), CD1C (cDC2), LAMP3 (mature) | liver | doi:10.1038/s41586-022-05400-x |
| Endothelial (incl. LSEC) | VWF, PLVAP; liver-sinusoidal CLEC4G, STAB2, OIT3, FCN3 | liver | doi:10.1038/s41586-022-05400-x |
| Fibroblast / HSC | COL1A1, ACTA2; pericyte/HSC RGS5 | liver | doi:10.1038/s41586-022-05400-x |
| Hepatocyte | ALB, APOA1, APOE, TTR, GPC3, AFP | liver parenchyma | doi:10.1038/s41586-022-05400-x |
| Cholangiocyte | EPCAM, KRT7, KRT19, SOX9, PROM1 | bile duct / ICC | doi:10.1038/s41586-022-05400-x |

**Liver-specific caveats (do not skip):**
- **Malignant vs normal epithelium:** in HCC, malignant hepatocytes retain ALB/HNF4A; in ICC,
  malignant cells are KRT19+/EPCAM+ — so `KRT19+` ≠ normal duct, and tumour cells can co-express
  hepatocyte + biliary markers. Separate malignant from normal with **inferred CNV (inferCNV)**,
  not markers alone. Expect strong patient-specific (clonal) separation of malignant epithelium.
- **Neutrophils / TANs** are fragile and low-mRNA — routinely lost in droplet scRNA-seq and nearly
  absent in snRNA-seq. S100A8/9 overlap heavily with classical monocytes; lean on FCGR3B / CSF3R /
  CXCR2 to separate, and don't over-interpret if few are recovered.
- Most liver endothelium is **sinusoidal (LSEC)** where VWF can be low — confirm with CLEC4G/STAB2.
