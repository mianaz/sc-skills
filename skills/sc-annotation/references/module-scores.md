# Module / signature scores (AddModuleScore)

Score a published gene signature per cell — TLS imprint, 12-chemokine, MHC-I/II, immunoglobulin,
etc. Use for state/program scoring on top of discrete cell-type labels.

```r
# AddModuleScore appends a column named <name><k>; rename to drop the trailing index.
seu <- Seurat::AddModuleScore(seu, features = list(tls_genes), name = "TLS")
seu$`TLS-imprint` <- seu$TLS1; seu$TLS1 <- NULL
```

## Build gene-family scores by pattern (mouse examples)
```r
mhc1 <- grep("^H2-(K|D|L|T|Q|M)", rownames(seu), value = TRUE)
mhc2 <- grep("^H2-(A|E)",          rownames(seu), value = TRUE)
ig   <- grep("^Ig[lkh]",           rownames(seu), value = TRUE)
for (s in list(c("MHC-I", "mhc1"), c("MHC-II", "mhc2"), c("Ig", "ig"))) {
  seu <- Seurat::AddModuleScore(seu, features = list(get(s[2])), name = s[1])
  seu[[s[1]]] <- seu[[paste0(s[1], "1")]]; seu[[paste0(s[1], "1")]] <- NULL
}
```

## Example published signatures (cite the source paper in methods.md)
- **TLS imprint** (Meylan et al., tertiary lymphoid structures): Ig genes (Igha/Ighg*/Igkc/Iglc*/Jchain),
  B markers (Cd79a, Mzb1, Xbp1, Fcrl5, Ssr4), T (Trbc2, Il7r), fibroblast (Cxcl12, Lum),
  complement (C1qa, C7), plus Cd52/Apoe/Pim2/Derl3.
- **12-chemokine TLS signature**: Ccl2/3/4/5/8/19/21, Cxcl9/10/11/13, Ccl17 (mouse orthologs;
  Ccl18 has no mouse ortholog — substitute Ccl8/Ccl17).

Scores are relative within a dataset. Compare across groups with a test + exact p-value
(sc-conventions); visualize as violin (tidyplots) or dotplot (scientific-plotting). For human signatures
in a mouse object, map orthologs first (sc-preprocessing/gene-symbols.md). UCell is an
alternative scorer that is rank-based and less sensitive to library size.

## Bridge: derive a signature from BULK, project onto single cells

A recurring cross-modality pattern (e.g. a bulk-defined cell-type/program signature scored per single cell). Direction bulk → sc:

```r
# 1. Signature from a bulk DE table (DESeq2/edgeR): top up-genes for the state/type of interest.
sig <- subset(bulk_de, padj < 0.05 & log2FoldChange > 1)
sig <- head(sig[order(-sig$log2FoldChange), "gene"], 100)      # ~50-200 genes
sig <- intersect(sig, rownames(seu))                            # map orthologs first if cross-species

# 2. Score per cell (UCell preferred here — rank-based, robust to the bulk/sc platform gap).
seu <- UCell::AddModuleScore_UCell(seu, features = list(bulk_sig = sig))
# then: violin/dotplot of bulk_sig_UCell across cell types; test enrichment in the expected type.
```

**Validate in the reverse direction** rather than trusting the projection alone: **pseudobulk the single-cell object** (→ `sc-pseudobulk`) and check that the signature genes' pseudobulk fold-changes agree with the bulk fold-changes (Spearman of logFCs, or a scatter). Concordant fold-changes + enrichment in the expected cell type = a trustworthy bridge; enrichment without fold-change concordance often means the signature is capturing library-size/ambient signal, not the biology. Cite the bulk source and record both checks in methods.md.
