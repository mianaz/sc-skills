# Volcano Plots

Use scop when the input is a Seurat/scop DE result or a DE table matching scop's expected columns. Use tidyplots or ggplot2 when the DE table is generic or needs custom layering.

## scop

```r
p <- scop::VolcanoPlot(
  srt = NULL,
  res = de_tbl,
  DE_threshold = "abs(avg_log2FC) > 1 & p_val_adj < 0.05",
  x_metric = "avg_log2FC",
  y_metric = "p_val_adj",
  features_label = c("ALB", "KRT19", "PECAM1"),
  nlabel = 8,
  palette = "RdBu"
)
```

Required table columns include `gene`, `group1`, `avg_log2FC`, and `p_val_adj`.

## tidyplots layered volcano

```r
library(dplyr)
library(tidyplots)

plot_df <- de_tbl |>
  mutate(
    gene = coalesce(gene, rownames(de_tbl)),
    neg_log10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
    direction = case_when(
      avg_log2FC >= 1 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC <= -1 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "NS"
    ),
    label_me = gene %in% genes_label
  )

p <- plot_df |>
  tidyplot(x = avg_log2FC, y = neg_log10_padj, width = 65, height = 55) |>
  add_data_points(data = filter_rows(direction == "NS"), color = "grey80", size = 0.5) |>
  add_data_points(data = filter_rows(direction == "Up"), color = "#B2182B", size = 0.7, alpha = 0.65) |>
  add_data_points(data = filter_rows(direction == "Down"), color = "#2166AC", size = 0.7, alpha = 0.65) |>
  add_reference_lines(x = c(-1, 1), y = -log10(0.05)) |>
  add_data_labels_repel(
    data = filter_rows(label_me),
    label = gene,
    background = TRUE,
    max.overlaps = Inf
  ) |>
  adjust_x_axis(title = "log2 fold change") |>
  adjust_y_axis(title = "-log10 adjusted p-value")
```

## ggplot2 fallback

```r
library(ggplot2)
library(cowplot)
library(ggrepel)

p <- ggplot(plot_df, aes(avg_log2FC, neg_log10_padj, color = direction)) +
  geom_point(size = 0.7, alpha = 0.65) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = subset(plot_df, label_me),
    aes(label = gene),
    box.padding = 0.35,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(NS = "grey80", Up = "#B2182B", Down = "#2166AC")) +
  theme_cowplot(font_size = 12)
```

Do not label every significant point. Label prespecified genes, top hits, or a small curated set.
