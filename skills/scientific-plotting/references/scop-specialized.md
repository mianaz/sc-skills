# scop Specialized Plot Families

Read this when the figure is not a basic UMAP/composition/marker plot. Many specialized plots expect that the matching `Run*()` function has already stored results in the Seurat object.

## Differential expression

```r
# From a Seurat object with stored DE, or directly from a result table.
p <- scop::VolcanoPlot(
  srt = NULL,
  res = de_tbl,
  DE_threshold = "abs(avg_log2FC) > 1 & p_val_adj < 0.05",
  features_label = c("ALB", "KRT19", "PECAM1"),
  nlabel = 10,
  y_metric = "p_val_adj",
  palette = "RdBu"
)

p2 <- scop::DEtestPlot(
  srt = NULL,
  res = de_tbl,
  plot_type = "manhattan",
  y_metric = "p_val_adj"
)
```

`VolcanoPlot()` accepts a DE table with columns including `gene`, `group1`, `avg_log2FC`, and `p_val_adj`. Use `features_label` for preselected genes; use `nlabel` sparingly.

## Enrichment, GSEA, GSVA, and pathway activity

```r
p_enrich <- scop::EnrichmentPlot(
  srt = seu,
  db = "GO_BP",
  group.by = "cell_type",
  plot_type = "dot",
  topTerm = 8,
  padjustCutoff = 0.05
)

p_gsea <- scop::GSEAPlot(
  srt = seu,
  db = "GO_BP",
  group.by = "cell_type",
  plot_type = "bar",
  direction = "both",
  topTerm = 8
)

ht_gsva <- scop::GSVAPlot(
  srt = seu,
  group.by = "cell_type",
  plot_type = "heatmap",
  topTerm = 30,
  show_row_names = TRUE
)
```

Use ggh4x or ggplot2 fallback when you need a custom pathway dotplot with wrapped names and colored source strips.

## Cell-cell communication

```r
p_dot <- scop::CCCHeatmap(
  seu,
  method = "CellChat",
  plot_type = "dot",
  sender.use = sender_cells,
  receiver.use = receiver_cells,
  top_n = 50
)

p_net <- scop::CCCNetworkPlot(
  seu,
  method = "CellChat",
  plot_type = "circle",
  signaling = "CXCL"
)

p_stat <- scop::CCCStatPlot(
  seu,
  method = "CellChat",
  plot_type = "bar",
  display_by = "interaction",
  top_n = 20
)
```

Use CCC plots only after the relevant communication workflow has run (`RunCCC()`, `RunCellChat()`, `RunCellphoneDB()`, `RunLIANA()`, etc.).

## Deconvolution, abundance, and differential abundance

```r
p_dec <- scop::DeconvolutionPlot(res = deconv_tbl, plot_type = "bar")

ht_dec <- scop::DeconvolutionPlot(
  res = deconv_tbl,
  plot_type = "heatmap",
  sample_annotation = "condition",
  sample_split = "condition"
)

p_immune <- scop::ImmuneAbundancePlot(
  object = seu,
  group.by = "condition",
  plot_type = "box",
  add_stat = TRUE
)

p_da <- scop::ProportionTestPlot(
  seu,
  plot_type = "effect",
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5)
)
```

Use bar/fill plots for composition; use boxplots for sample-level abundance comparisons; use effect plots for differential abundance results.

## Spatial gradients and spot composition

```r
p_grad <- scop::SpatialGradientPlot(
  seu,
  plot_type = "combined",
  features = c("COL1A1", "MKI67"),
  overlay_image = TRUE
)

p_spot_pie <- scop::SpatialSpotPlot(
  seu,
  plot_type = "pie",
  group.by = "RCTD_dominant_type",
  overlay_image = TRUE
)
```

Use this reference for final rendering of spatial plots.

### Visium spatial dual-layer plot (continuous fill + categorical overlay)
(source: Cho et al. Science 2026)

For plots that show a continuous score across all spots AND highlight selected spots with
a categorical marker (e.g. colocalization heatmap + candidate TLS spots). Uses ggplot2
directly, not scop.

```r
library(ggplot2); library(RColorBrewer)
pal <- colorRampPalette(brewer.pal(9, "GnBu"))

ggplot(all_spots, aes(x = col, y = -row, fill = score)) +
  geom_point(shape = 21, colour = "white", size = 1.5, stroke = 0.5) +
  geom_point(data = selected_spots, aes(x = col, y = -row),
             pch = 5, fill = NA, size = 1.5, colour = "#00afb9", stroke = 0.3) +
  scale_fill_gradientn(colours = pal(100), limits = c(0, max_score)) +
  labs(x = NULL, y = NULL, fill = "Score") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_blank(), axis.ticks = element_blank())
```

Key choices: shape=21 (filled circle w/ border) for background, pch=5 (diamond) for
overlay to visually separate layers; flip Y (`y = -row`) to match tissue orientation;
clean theme with no gridlines preserves spatial context.

## Metabolism, TF activity, GRN, and ATAC tracks

```r
p_tf <- scop::DorotheaPlot(
  seu,
  group.by = "condition",
  group1 = "treated",
  group2 = "control",
  top_n = 25,
  p.adjust.method = "BH"
)

p_met <- scop::MetabolismPlot(
  seu,
  group.by = "cell_type",
  plot_type = "heatmap",
  topTerm = 20,
  show_row_names = TRUE
)

p_scenic <- scop::SCENICPlot(
  seu,
  group.by = "cell_type",
  plot_type = "rss_dotplot",
  top_n = 12
)

p_atac <- scop::CoverageTrackPlot(
  seu,
  region = "chr1-1000000-1010000",
  assay = "peaks",
  group.by = "cell_type",
  links = TRUE
)
```

When the specialized plot returns a ComplexHeatmap object or list, explicitly draw/save the plot object rather than assuming `ggsave()` will work.
