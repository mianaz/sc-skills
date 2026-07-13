# Curated reference gene sets

Reusable gene sets for AddModuleScore, ssGSEA, or marker-based detection. Each carries
provenance so the source paper can be checked.

## TLS (tertiary lymphoid structure) marker genes
(source: Cho et al. Science 2026)

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
with spatial colocalization for TLS boundary detection.
