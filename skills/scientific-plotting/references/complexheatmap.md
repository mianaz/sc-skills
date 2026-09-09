# ComplexHeatmap — annotated heatmaps (row/col annotations, bar annotations, column split)

Use when pheatmap isn't enough: side color-bar annotations, an annotation that is itself a bar plot,
or splitting columns/rows into labeled blocks. House RdBu palette + white tile borders.

```r
library(ComplexHeatmap); library(circlize); library(RColorBrewer)

colors    <- rev(brewer.pal(11, "RdBu"))
colors.use<- grDevices::colorRampPalette(colors)(100)
my_breaks <- c(seq(-2, 0, length.out = 50), seq(0.01, 2, length.out = 50))
col_fun   <- colorRamp2(breaks = my_breaks, colors = colors.use)   # continuous mapping

# LEFT annotation: categorical group color bar (rows = samples here)
row_ha <- rowAnnotation(
  Group = sample_groups$group,
  col = list(Group = c(DKO = "#E41A1C", TKO = "#377EB8", TKO_CSI = "#4DAF4A")),
  width = unit(5, "mm"), show_annotation_name = FALSE)

# RIGHT annotation: a bar plot of a per-row numeric value
ratio_bar <- rowAnnotation(
  "CD8/Treg\nratio" = anno_barplot(cell_counts$cd8_treg_ratio, bar_width = 0.8,
                                   gp = gpar(fill = "#7B1FA2", col = "white"),
                                   width = unit(2.5, "cm")),
  annotation_name_rot = 0)

# column blocks (e.g. CD8 vs CD4 gene sets)
col_block <- factor(c(rep("CD8 T cells", n_cd8), rep("CD4 T cells", n_cd4)),
                    levels = c("CD8 T cells", "CD4 T cells"))

ht <- Heatmap(mat, name = "Z-score", col = col_fun,
              cluster_rows = FALSE, cluster_columns = FALSE,
              column_split = col_block, column_gap = unit(3, "mm"),
              rect_gp = gpar(col = "white", lwd = 0.5), border = TRUE,
              row_names_gp = gpar(fontsize = 14), column_names_gp = gpar(fontsize = 14))

ht_list <- row_ha + ht + ratio_bar     # compose annotations around the heatmap

png("figures/heat_clean.png", width = 10, height = 3.2, units = "in", res = 600)
draw(ht_list); dev.off()
pdf("figures/heat_clean.pdf", width = 10, height = 3.2); draw(ht_list); dev.off()
```

`column_split` (and `row_split`) accept a factor to make labeled blocks; `anno_barplot`/`anno_points`/
`anno_boxplot` build data-driven annotations; combine annotations with `+` (rows) or stacking
(`%v%`, columns). Use ComplexHeatmap when you need these; otherwise the simpler decoupleR pheatmap
(references/heatmaps.md) is fine. Save dual-version + the matrix as source data (sc-conventions).

## Row-split with category annotation sidebar
(source: sc-paper-distill/papers/2026-pan-cancer-TLS.md F2)

For marker-by-cluster heatmaps where genes belong to functional categories (e.g. T cells /
B cells / Initiating / Activated), split rows by category and add a color-coded sidebar.

```r
library(ComplexHeatmap); library(circlize)

# 3-point color ramp (blue → yellow → red)
col_fun <- colorRamp2(c(-1, 0, 5), c("#118ab2", "#fdffb6", "#e63946"))

# Row category annotation
category_colors <- c("T cells" = "#e76f51", "B cells" = "#00CDAC",
                     "Initiating" = "#f4a261", "Activated" = "#0077b6")
row_cat <- factor(c(rep("T cells", 2), rep("B cells", 4),
                    rep("Initiating", 2), rep("Activated", 26)),
                  levels = names(category_colors))
ha <- HeatmapAnnotation(
  df = data.frame(Category = row_cat), which = "row",
  col = list(Category = category_colors))

Heatmap(zscored_mat, col = col_fun,
        cluster_columns = FALSE, cluster_rows = FALSE,
        show_column_names = FALSE, show_row_names = TRUE,
        column_names_side = "top", row_names_side = "right",
        row_split = row_cat, column_split = col_split_factor,
        row_gap = unit(3, "mm"), column_gap = unit(3, "mm"),
        left_annotation = ha,
        heatmap_legend_param = list(title = "Z-score",
                                    legend_height = unit(3, "cm")))
```
