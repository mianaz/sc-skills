# scop trajectory / method-specific plots — run THEN plot

Method-specific scop plots require running scop's built-in analysis function first
so the results live in the object in the form the plot expects.

| Analysis (run first) | Plot function(s) | Notes |
|---|---|---|
| RunSlingshot | LineagePlot; CellDimPlot(lineages=) | LineagePlot draws lineage curves; or overlay on CellDimPlot via the `lineages` arg |
| RunPAGA | PAGAPlot | PAGA result stored in `srt@misc$paga`; type = connectivities / connectivities_tree |
| RunDynamicFeatures | DynamicHeatmap; DynamicPlot | pseudotime heatmap of dynamic features along lineages |
| RunSCVELO | VelocityPlot | RNA velocity; plot_type = raw / grid / stream |
| RunCytoTRACE | CytoTRACEPlot | differentiation potential / pseudotime potency |
| RunPalantir / RunCellRank | PseudotimeProjectionPlot; BranchStreamPlot | pseudotime / fate projections |

## Slingshot — run, then plot
```r
seu <- scop::RunSlingshot(seu, group.by = "cell_type", reduction = "umap")
# dedicated lineage-curve plot:
p   <- scop::LineagePlot(seu, lineages = c("Lineage1","Lineage2"), reduction = "umap")
# or overlay lineages on a dim plot:
p   <- scop::CellDimPlot(seu, group.by = "cell_type", reduction = "umap",
                         lineages = c("Lineage1","Lineage2"))
```

## PAGA — run, then plot
```r
seu <- scop::RunPAGA(seu, group.by = "cell_type", reduction = "umap")
p   <- scop::PAGAPlot(seu, reduction = "umap")
```

## Dynamic features (pseudotime heatmap) — run, then plot
```r
seu <- scop::RunSlingshot(seu, group.by = "cell_type", reduction = "umap")
seu <- scop::RunDynamicFeatures(seu, lineages = c("Lineage1","Lineage2"))
p   <- scop::DynamicHeatmap(seu, lineages = c("Lineage1","Lineage2"))
```

## RNA velocity — run, then plot
```r
seu <- scop::RunSCVELO(seu, group.by = "cell_type", reduction = "umap")
p   <- scop::VelocityPlot(seu, reduction = "umap", plot_type = "stream")
```
