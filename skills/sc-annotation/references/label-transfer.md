# Label transfer from a curated reference atlas
```r
anchors <- Seurat::FindTransferAnchors(reference = ref, query = seu,
                                       dims = 1:30, reference.reduction = "pca")
pred <- Seurat::TransferData(anchorset = anchors, refdata = ref$cell_type, dims = 1:30)
seu$transfer_label <- pred$predicted.id
seu$transfer_score <- pred$prediction.score.max
```
Always cite the reference atlas (PMID/DOI) in methods.md. Treat transferred labels as CANDIDATES — confirm with markers (references/marker-verification.md).
