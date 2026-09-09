# Heatmaps

Use the simplest heatmap that answers the question. Save the exact plotted matrix as source data.

## pheatmap for compact matrices

Use this for pathway/activity matrices, marker z-scores, and simple clustered heatmaps.

```r
colors <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
colors_use <- grDevices::colorRampPalette(colors)(100)
breaks <- c(seq(-2, 0, length.out = 51), seq(0.04, 2, length.out = 50))

readr::write_csv(as.data.frame(mat) |> tibble::rownames_to_column("feature"),
                 "source_data/heatmap_matrix.csv")

png("figures/heatmap_clean.png", width = 7, height = 5, units = "in", res = 600)
pheatmap::pheatmap(
  mat,
  color = colors_use,
  breaks = breaks,
  border_color = "white",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45
)
dev.off()

pdf("figures/heatmap_clean.pdf", width = 7, height = 5)
pheatmap::pheatmap(mat, color = colors_use, breaks = breaks, border_color = "white")
dev.off()
```

## tidyplots heatmap for tidy tables

```r
p <- heat_df |>
  tidyplots::tidyplot(x = condition, y = feature, color = value, width = 55, height = 70) |>
  tidyplots::add_heatmap(scale = "row", rotate_labels = 45, rasterize = TRUE) |>
  tidyplots::adjust_colors(tidyplots::colors_diverging_BuRd)
```

## Diagonalize gene order (readable marker matrices)

For a cell-type × marker matrix, input gene order is rarely readable and clustering genes
scatters each type's markers. **Diagonalize** instead: order genes so high expression runs
down the diagonal. Two strategies — pick by data shape.
(source: sc-paper-distill/papers/2019-FCA-liver.md F3 — Popescu et al. Nature 2019; papers/2020-descartes-fetal-atlas.md F3, F4 — Cao et al. Science 2020, argmax staircase + Z-cap.)

**Scale before diagonalizing: Z-score each gene across types, then cap** (Cao capped to `[0, 3]`).
Without a cap, one high-expression gene sets the color scale and flattens everything else to a single
hue; capping the Z-score keeps the block-diagonal signature legible. The same capped Z-scale should
drive the companion dotplot's color so heatmap and dotplot read identically.

```r
# mat: rows = cell types (in your biological/lineage order), cols = genes; values = mean expr.

# (A) Weighted center-of-mass — smooth diagonal, best for graded expression.
weighted_center <- function(v) sum(seq_along(v) * v) / sum(v)
mat <- mat[, order(apply(mat, 2, weighted_center))]

# (B) Argmax + block — sharp blocks, best for discrete type-specific markers.
z <- scale(mat); z[is.na(z)] <- 0          # z-score each gene across cell types
owner <- apply(z, 2, which.max)            # the cell type each gene peaks in
mat   <- mat[, order(owner, -apply(z, 2, max))]
# owner_block <- factor(rownames(mat)[owner][order(owner, ...)])  # -> column_split in ComplexHeatmap

# (C) Alternative: cluster genes by similarity (not diagonalized)
# mat <- mat[, hclust(dist(t(mat)), method = "ward.D")$order]
```

Reuse the **same** diagonalized order across the heatmap and its companion dotplot/spotplot so
reviewers can map one to the other. Save the computed order with the source data so the figure
is reproducible. For the block version, feed `owner_block` to ComplexHeatmap `column_split`.

## pheatmap with significance stars overlay
(source: sc-paper-distill/papers/2026-pan-cancer-TLS.md F3)

For score-by-condition matrices with statistical testing (e.g. ICB response across
cohorts, pathway activity across conditions), overlay significance stars on cells.

```r
library(pheatmap)

# Custom diverging palette (dark red → white → dark blue)
pal <- rev(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7',
             '#ffe8d6','#D1E5F0','#4393C3','#2166AC','#053061'))

# Build significance matrix from p-values
pmt <- matrix(nrow = nrow(p_matrix), ncol = ncol(p_matrix),
              dimnames = dimnames(p_matrix))
pmt[p_matrix < 0.01]                      <- '**'
pmt[p_matrix >= 0.01 & p_matrix < 0.05]   <- '*'
pmt[p_matrix >= 0.05]                      <- ''

pheatmap(score_matrix, scale = "none",
         cluster_rows = TRUE, cluster_cols = TRUE, border = NA,
         display_numbers = pmt, fontsize_number = 12, number_color = "white",
         cellwidth = 10, cellheight = 10, color = pal)
```

## Cross-validation confusion matrix (label/cluster validity)
To show that annotations are separable (not over-clustered), plot the **predicted × actual recall
matrix** from a cross-validated classifier as a heatmap, ordered so the diagonal reads, paired with an
F1 box plot vs a permuted-label null. High on-diagonal recall + F1 ≫ permuted = trustworthy labels.
(source: sc-paper-distill/papers/2020-descartes-fetal-atlas.md F2 — Cao et al. Science 2020; classifier recipe in sc-annotation/references/cluster-specificity.md)

```r
# cm: rows = predicted, cols = actual (row-normalized to recall, 0..1). Order rows/cols the same way.
ord <- hclust(dist(cm))$order          # or a fixed biological order shared with other panels
cm  <- cm[ord, ord]
readr::write_csv(as.data.frame(cm) |> tibble::rownames_to_column("predicted"),
                 "source_data/cv_confusion.csv")

pheatmap::pheatmap(cm, cluster_rows = FALSE, cluster_cols = FALSE,
                   color = viridisLite::viridis(100), border_color = "white",
                   cellwidth = 10, cellheight = 10,        # square tiles
                   display_numbers = FALSE, angle_col = 45)
```

- **Viridis 0→1, square tiles, no clustering** once ordered — the message is the diagonal, so fix the
  order (shared with the companion dotplot/UMAP) rather than re-clustering per panel.
- **Recall (row-normalized), not raw counts:** counts let the biggest class dominate the color scale;
  recall makes every type's self-recovery comparable. Off-diagonal blocks name the confusable pairs.
- Companion **F1 box plot** (permuted null vs main types vs subtypes) belongs in tidyplots — see
  `references/tidyplots.md`; always include the permuted-label box as the floor.

## Escalate to ComplexHeatmap when needed

Use `references/complexheatmap.md` when you need:

- row/column annotations
- split blocks with labels
- barplot/point annotations beside the heatmap
- separate legends
- explicit row/column name control
- drawing/saving a `Heatmap` object
