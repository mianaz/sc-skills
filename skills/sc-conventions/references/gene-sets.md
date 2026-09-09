# Curated reference gene sets

Reusable gene sets for AddModuleScore, ssGSEA, or marker-based detection. Each carries
provenance so the source paper can be checked.

## TLS (tertiary lymphoid structure) marker genes
(source: sc-paper-distill/papers/2026-pan-cancer-TLS.md §4 — Cho et al. Science 2026)

For detecting TLS regions in any tissue via Seurat `AddModuleScore` on spatial or
scRNA-seq data. These genes mark B-cell-rich lymphoid aggregates associated with
antitumor immunity.

```r
tls_markers <- c("MS4A1", "CXCR5", "SELL", "CD19", "LTB",
                 "CD79B", "CD37", "CD79A", "TCL1A")

seu <- AddModuleScore(seu, features = list(
  intersect(tls_markers, rownames(seu))), nbin = 10, name = "TLS_score")
```

Applicable to: any human tissue with potential immune aggregates (liver, lung, breast,
pancreas, etc.). The score identifies B-cell-dominated lymphoid structures; combine
with spatial colocalization (sc-spatial) for TLS boundary detection.

## Courau conserved TME gene expression programs
(source: sc-paper-distill/papers/2026-humu-tme.md M5–M6 — Courau et al. Nat Immunol 2026)

Cross-species conserved GEPs. Score T_3 in T cells and My_2 in myeloid cells, then test their
coordination (a “movement”) and 2×2 survival. Full ranked lists are in the paper’s Supplementary
Table 2; genes below are the drivers named in the main text.

```r
courau_t3_cytotoxicity <- c("PRF1", "LAG3", "GZMB", "NKG7")
courau_t9_cd4_regulation <- c("TSC22D3", "JUNB", "RGS1", "IL7R", "CD69")
courau_my1_inflammatory <- c("IL1A", "IL1B", "NLRP3")
courau_my2_ifn <- c("IFIT2", "IFIT3", "ISG15", "CXCL10")
courau_my11_lyve1_tam <- c("LYVE1", "FOLR2", "CD163", "MRC1", "SELENOP")
```

Do not use CXCL13+ T / TLS biology in standard mouse models without an explicit species caveat
(human CXCL13 is T/Treg-biased; mouse CXCL13 is typically stromal).
