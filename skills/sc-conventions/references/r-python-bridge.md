# R ↔ Python bridge

Seurat is the source of truth. Export to AnnData/h5ad for Python steps, then
carry results back. Preserve **raw counts in a layer** (scVI needs counts).

## Export Seurat → h5ad
Preferred: scConvert (github: mianaz/scConvert) or anndataR; SeuratDisk as fallback.
```r
# SeuratDisk fallback (widely available)
library(SeuratDisk)
SaveH5Seurat(seu, filename = "objects/seu.h5Seurat", overwrite = TRUE)
Convert("objects/seu.h5Seurat", dest = "h5ad", overwrite = TRUE)
# -> objects/seu.h5ad ; ensure counts are in a layer for scVI
```
```r
# scConvert (mianaz/scConvert): universal converter, dest sets output format.
scConvert(seu, dest = "objects/seu.h5ad")
# reverse direction: scConvert("objects/seu.h5ad", dest = "h5seurat")

# anndataR (scverse/anndataR): convert Seurat -> AnnData, then write h5ad.
adata <- anndataR::as_AnnData(seu)
anndataR::write_h5ad(adata, "objects/seu.h5ad")
# read back with conversion: anndataR::read_h5ad("objects/seu.h5ad", as = "Seurat")
```

## Carry embeddings back into Seurat
```r
# latent: cells x dims matrix exported from Python (e.g. read.csv of obsm['X_scVI'])
rownames(latent) <- colnames(seu)            # align to Seurat cell order
colnames(latent) <- paste0("scVI_", seq_len(ncol(latent)))
seu[["scvi"]] <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(latent), key = "scVI_", assay = DefaultAssay(seu))
```

## Carry labels/scores back
```r
seu$scanvi_pred <- preds[colnames(seu)]       # align by cell barcode
```
Always align by cell barcode, never by position assumption.
