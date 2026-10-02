# Per-cluster regulon summary (binarize + activity heatmap)

After AUCell (pipeline stage 3), to report which regulons mark which clusters, binarize the AUCell matrix and aggregate per cluster.

(source: https://doi.org/10.1038/s41586-020-2316-7)

```r
# bin: binary regulon-activity matrix (regulons × cells), from AUCell_exploreThresholds()/binarize.
cor_m <- cor(t(bin), method = "spearman")
bin   <- bin[rowSums(cor_m > 0.3) <= 1, ]                 # drop redundant regulons
act   <- sapply(split(colnames(bin), clusters),           # per-cluster mean activity
                \(cl) rowMeans(bin[, cl, drop = FALSE]))
act   <- act[rowSums(act >= 0.3) >= 1, ]                  # keep regulons active in ≥30% of ≥1 cluster
pheatmap::pheatmap(act, clustering_method = "ward.D2", border_color = NA)
```

Binarization (GMM/bimodal thresholding) turns AUCell's relative scores into on/off calls that compare better across clusters. decoupleR (`references/decoupler-activity.md`) is the motif-DB-free alternative.
