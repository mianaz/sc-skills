# SingleR + celldex (reference-based candidate labels)

Correlation-based automated labelling against a curated bulk/sc reference. Fast, no anchor
finding. Good complement to label transfer and LLM calls. **Candidates only** — verify with markers.

```r
library(SingleR); library(celldex)

ref <- celldex::ImmGenData()            # mouse immune; others: MouseRNAseqData,
                                        # HumanPrimaryCellAtlasData, BlueprintEncodeData,
                                        # DatabaseImmuneCellExpressionData, MonacoImmuneData
pred <- SingleR(test    = seu[["RNA"]]$data,   # log-normalized
                ref     = ref,
                labels  = ref$label.fine,      # or ref$label.main for coarse
                clusters = seu$seurat_clusters)# cluster-level (omit for per-cell)

# transfer cluster-level labels back to cells
lab <- setNames(pred$labels, rownames(pred))
seu$celltype_singler <- lab[as.character(seu$seurat_clusters)]
```

Pick the reference matching tissue/species; cite it (PMID/DOI) in methods.md. Cluster-level
labelling (`clusters=`) is more stable than per-cell. Treat as candidates; the canonical-marker
dotplot remains the mandatory gate (marker-verification.md).
