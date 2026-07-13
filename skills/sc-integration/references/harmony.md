# Harmony (R, fast path)
```r
seu <- Seurat::NormalizeData(seu) |>
       Seurat::FindVariableFeatures() |>
       Seurat::ScaleData() |>
       Seurat::RunPCA(npcs = 30)
seu <- harmony::RunHarmony(seu, group.by.vars = "sample")   # batch/sample key; saves reduction "harmony"
seu <- Seurat::RunUMAP(seu, reduction = "harmony", dims = 1:30)
seu <- Seurat::FindNeighbors(seu, reduction = "harmony", dims = 1:30) |>
       Seurat::FindClusters()
```
Use for fast exploration. Log the batch key + dims in methods.md.
