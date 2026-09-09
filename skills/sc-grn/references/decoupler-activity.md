# decoupleR — fast R-native TF & pathway activity (no motif DBs, no ATAC)

When you want per-cell TF or pathway activity in minutes from an existing Seurat object — and
do NOT want GRNBoost2/cisTarget/AUCell, 20GB motif databases, or scATAC — use decoupleR with a
curated prior network. This scores activity by enrichment of a TF/pathway's known targets in each
cell, not de-novo network inference.

## TF activity — CollecTRI + univariate linear model (run_ulm)
```r
library(decoupleR)
net <- get_collectri(organism = "mouse", split_complexes = FALSE)   # curated TF->target
mat <- seu[["SCT"]]$data |> as.matrix()                              # genes x cells, normalized
acts <- run_ulm(mat = mat, net = net,
                .source = "source", .target = "target", .mor = "mor", minsize = 5)

# store as an assay so you can DimPlot/FeaturePlot TF activity like genes
seu[["tfsulm"]] <- acts |>
  tidyr::pivot_wider(id_cols = "source", names_from = "condition", values_from = "score") |>
  tibble::column_to_rownames("source") |>
  Seurat::CreateAssayObject()
DefaultAssay(seu) <- "tfsulm"
seu <- Seurat::ScaleData(seu)
seu[["tfsulm"]]$data <- seu[["tfsulm"]]$scale.data
```

## Pathway activity — PROGENy + multivariate linear model (run_mlm)
```r
net <- get_progeny(organism = "mouse", top = 1000)                   # 14 signaling pathways
acts <- run_mlm(mat = mat, net = net,
                .source = "source", .target = "target", .mor = "weight", minsize = 5)
seu[["pathwaysmlm"]] <- acts |>
  tidyr::pivot_wider(id_cols = "source", names_from = "condition", values_from = "score") |>
  tibble::column_to_rownames("source") |>
  Seurat::CreateAssayObject()
# (scale as above)
```

Use `organism = "human"` for human data. `run_ulm` for the large sparse CollecTRI net; `run_mlm`
suits the dense weighted PROGENy net. Scores are relative within a dataset (compare cells/groups
here, not across datasets). Write the long `acts` table to a CSV and log the prior + method in methods.md.

## Mean-activity heatmap per group (cell-type × condition)
Aggregate the activity assay to group means, pick the most variable sources, plot RdBu square tiles
(sc-conventions palette / scientific-plotting decoupleR heatmap):
```r
df <- t(as.matrix(seu@assays$tfsulm$data)) |> as.data.frame() |>
  dplyr::mutate(group = paste0(seu$celltype, "_", seu$condition)) |>
  tidyr::pivot_longer(-group, names_to = "source", values_to = "score") |>
  dplyr::group_by(group, source) |> dplyr::summarise(mean = mean(score), .groups = "drop")
top <- df |> dplyr::group_by(source) |> dplyr::summarise(sd = sd(mean)) |>
  dplyr::slice_max(sd, n = 25) |> dplyr::pull(source)
mat_g <- df |> dplyr::filter(source %in% top) |>
  tidyr::pivot_wider(names_from = source, values_from = mean) |>
  tibble::column_to_rownames("group") |> as.matrix()
# pheatmap with rev(brewer.pal(11,"RdBu")), symmetric breaks centered at 0 (scientific-plotting/heatmaps.md)
```
