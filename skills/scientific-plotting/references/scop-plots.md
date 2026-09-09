# scop Core Plots

Use `scop` for Seurat-centric single-cell, spatial, and omics plots. This is `mengxu98/scop`, not the unmaintained `SCP`. Installed local version observed: `scop` 0.8.9. Always check `?scop::FunctionName` when using an unfamiliar argument.

`scop` usually takes a Seurat object as `srt`. Many plot functions return ggplot/patchwork objects; heatmap-like functions can return ComplexHeatmap objects.

## Core function map

| Need | Function |
|---|---|
| UMAP/tSNE/PCA groups | `CellDimPlot()` |
| 3D embedding | `CellDimPlot3D()` |
| Feature on embedding | `FeatureDimPlot()`, `FeatureDimPlot3D()` |
| Cell density on embedding | `CellDensityPlot()` |
| Composition and metadata stats | `CellStatPlot()` |
| Marker dot/violin/box/col/bar | `FeatureStatPlot()` |
| Marker/group expression heatmap | `GroupHeatmap()`, `FeatureHeatmap()` |
| Feature or cell correlation | `FeatureCorPlot()`, `CellCorHeatmap()` |
| Query/reference projection | `ProjectionPlot()` |
| Spatial spots and spatial gradients | `SpatialSpotPlot()`, `SpatialGradientPlot()` |
| Integration diagnostics | `DimsEstimatePlot()`, `LISIPlot()`, `BenchmarkPlot()` |

## Reusing SCP / scplotter code with the scop fork

Lots of published figure code uses **SCP** or **scplotter** (same plot family, e.g. `CellDimPlot`,
`FeatureDimPlot`). Our fork is `mengxu98/scop`; the *style values transfer*, only some arg names
differ. (source: sc-paper-distill/papers/2026-human-pregastrula.md F1, M1 — Nature 2026)

| SCP / scplotter | scop fork |
|---|---|
| `group_by` | `group.by` |
| `split_by` | `split.by` |
| `theme = "theme_blank"` | `theme_use = "theme_blank"` |
| `palcolor`, `pt.size`, `raster`, `label`, `legend.position`, `reduction` | unchanged |

```r
# published (SCP/scplotter):  CellDimPlot(obj, group_by="celltype", theme="theme_blank", split_by="sample", ...)
# scop fork:
scop::CellDimPlot(srt, group.by = "celltype", reduction = "umap", label = FALSE,
                  palcolor = mycols, pt.size = 0.01, theme_use = "theme_blank",
                  legend.position = "right", raster = FALSE, split.by = "sample")
```
Dense-UMAP defaults from that paper: `theme_use="theme_blank"`, `pt.size=0.01`, `raster=FALSE`
(see sc-conventions/references/figure-sizing-and-type.md for pt.size-by-density + save dims).
scop's first arg is `srt`; verify any unfamiliar arg with `?scop::CellDimPlot`.

## UMAP or embedding by group

```r
library(scop)

p <- scop::CellDimPlot(
  seu,
  group.by = "cell_type",
  reduction = "umap",
  label = length(unique(seu$cell_type)) <= 12,
  theme_use = "theme_blank",
  pt.size = 0.2
)
```

If labels collide, set `label = FALSE`, split by condition, or label only a coarser metadata column. Do not keep overlapping cluster names.

## Feature expression on embedding

```r
p <- scop::FeatureDimPlot(
  seu,
  features = c("ALB", "KRT19", "PECAM1"),
  reduction = "umap",
  theme_use = "theme_blank",
  combine = TRUE,
  ncol = 3
)
```

Use the same feature limits across panels when making a comparative figure.

For an annotation **marker panel**, group the features by lineage and make one panel per
compartment (T/NK, B/plasma, myeloid, stromal, parenchymal) rather than one giant grid —
it doubles as the annotation rationale. Curated liver/TIME marker sets live in
`sc-annotation/references/markers.md`. Seurat-native fallback when scop is unavailable:
`Seurat::FeaturePlot(seu, features, ncol = 3, label = TRUE)`.
(source: sc-paper-distill/papers/2022-scPLC.md F2)

For a top-N marker dotplot straight from `FindAllMarkers`, prefer `FeatureStatPlot(plot_type="dot")`
below; the Seurat-native fallback is `DotPlot(seu, features) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))`.
NOTE: Seurat ≥4 uses `avg_log2FC` — older code (and some published scripts) use `avg_logFC`. (F3)

## One embedding, recolored N ways (+ grey focal overlay, + trajectory arrows)
For atlas-scale figures, compute **one** UMAP and recolor the *same coordinates* by several metadata
columns (species, stage, tissue, cell type) as a small-multiple — the reader compares panels directly
because the geometry is fixed. Grey out the context and color only the focal set to spotlight one group.
(source: sc-paper-distill/papers/2020-descartes-fetal-atlas.md F6, F5 — Cao et al. Science 2020)

```r
# same embedding recolored by several columns (identical coords across panels)
keys <- c("species", "stage", "tissue", "cell_type")
ps <- lapply(keys, function(k) scop::CellDimPlot(seu, group.by = k, reduction = "umap",
                                                 theme_use = "theme_blank", pt.size = 0.1))
patchwork::wrap_plots(ps, ncol = 2)

# grey background = context, color = focal set (e.g. one species / one organ)
seu$focal <- ifelse(seu$species == "human", as.character(seu$cell_type), NA)   # NA -> grey
scop::CellDimPlot(seu, group.by = "focal", reduction = "umap",
                  na_color = "grey90", theme_use = "theme_blank", pt.size = 0.1)
```

- **Why fixed coords:** recoloring one embedding (vs. re-running UMAP per subset) keeps positions
  comparable across panels — differences you see are label differences, not layout differences.
- **Directional trajectory arrows** are a hand-annotation overlay on the UMAP (they encode inferred flow,
  e.g. HSPC → lineages); add with `annotate("segment", ..., arrow = arrow())` in a ggplot layer, and say
  in the caption they are drawn from the trajectory model, not computed vectors.
- Keep one shared palette/legend across the recolored panels (export it once — see
  `scientific-reproducibility/references/figure-output-contract.md`, separate-legend).

## Labeled per-group UMAP small-multiple grid (atlas "table of contents")
A grid of per-organ (or per-compartment) UMAPs with **on-plot cluster text labels and no side legend** —
each panel is self-documenting, and an empty grid cell can hold summary stats.
(source: sc-paper-distill/papers/2020-descartes-fetal-atlas.md F1 — Cao et al. Science 2020)

```r
organs <- levels(factor(seu$organ))
panels <- lapply(organs, function(o) {
  s <- seu[, seu$organ == o]
  scop::CellDimPlot(s, group.by = "cell_type", reduction = "umap",
                    label = TRUE, label_insitu = TRUE,   # text on clusters, not a legend
                    legend.position = "none", theme_use = "theme_blank", pt.size = 0.1) +
    ggplot2::ggtitle(o)
})
patchwork::wrap_plots(panels, ncol = 4)   # fill a spare cell with a summary-stats textGrob
```

- **Why in-plot labels, no legend:** at 15+ panels a shared legend is unreadable and a per-panel legend
  wastes space; labeling clusters in situ makes each panel stand alone. If labels collide, coarsen the
  labeled column or label only the largest types (see layout-legibility.md).
- Verify `label_insitu`/label args against `?scop::CellDimPlot` (scop fork arg names differ from SCP).

## Composition and metadata stats

```r
p_pct <- scop::CellStatPlot(
  seu,
  stat.by = "cell_type",
  group.by = "condition",
  plot_type = "bar",
  stat_type = "percent",
  position = "stack"
)

p_count <- scop::CellStatPlot(
  seu,
  stat.by = "cell_type",
  group.by = "condition",
  plot_type = "bar",
  stat_type = "count",
  position = "stack"
)
```

`CellStatPlot()` also supports plot styles such as bar, rose, ring, pie, trend, area, dot, and sankey. Use percent for composition claims; use count when sample size/cell yield is the point.

## Marker verification dotplot (genes × cell types)

**Use `GroupHeatmap(add_dot = TRUE)`** for canonical-marker panels.  
Do **not** use `FeatureStatPlot(plot_type = "dot")` — that facets **one gene per panel**.

```r
markers <- c("ALB", "APOA1", "KRT19", "PECAM1", "PTPRC")

ht <- scop::GroupHeatmap(
  seu,
  features = markers,
  group.by = "cell_type",
  add_dot = TRUE,
  flip = TRUE,                 # genes on x, cell types on y
  nlabel = 0,                  # avoid duplicate anno_mark labels
  row_title = character(0),    # suppress group.by name as y title when flip=TRUE
  column_names_rot = 45,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  layer = "data",
  exp_method = "zscore",
  lib_normalize = FALSE
)
# Save ht$plot with a device larger than heatmap width/height + plot.margin
```

## Violins / boxes / bars (per-feature by group)

```r
p_vln <- scop::FeatureStatPlot(
  seu,
  stat.by = markers,
  group.by = "cell_type",
  plot_type = "violin"
)
```

Use tidyplots for non-Seurat tidy expression tables.

## Group heatmaps (no dots)

```r
p <- scop::GroupHeatmap(
  seu,
  features = markers,
  group.by = "cell_type"
)
```

For annotated matrices with custom side bars or split blocks, switch to `ComplexHeatmap`.

## Density and projections

```r
p_density <- scop::CellDensityPlot(
  seu,
  features = "module_score",
  group.by = "condition",
  split.by = "condition"
)

p_proj <- scop::ProjectionPlot(
  srt_query = query,
  srt_ref = reference,
  query_group = "predicted_cell_type",
  ref_group = "cell_type"
)
```

Use density plots when overplotting hides whether a group occupies a region. Use projection plots for reference mapping or label-transfer diagnostics.

## Spatial spots

```r
p_spatial <- scop::SpatialSpotPlot(
  seu,
  group.by = "cell_type",
  image = NULL,
  overlay_image = TRUE,
  pt.alpha = 0.9
)

p_gene <- scop::SpatialSpotPlot(
  seu,
  features = "MKI67",
  overlay_image = TRUE,
  palette = "Spectral"
)
```

Use `plot_type = "pie"` for spot-level cell-type proportions when deconvolution output is stored as numeric metadata columns.

## Integration diagnostics

```r
p_lisi <- scop::LISIPlot(
  seu,
  features = c("batch_LISI", "celltype_LISI"),
  reduction = "umap",
  plot_boxplot = TRUE
)
```

Use diagnostics to show batch mixing and biological conservation after integration.
