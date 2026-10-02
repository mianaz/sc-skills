# R-native pseudotime (Seurat-first, no Python)

When you work in R and don't want the scanpy route, two lightweight options:

(source: https://doi.org/10.1038/s41586-020-2316-7)

```r
# Diffusion pseudotime in R (destiny) — equivalent to scanpy dpt.
pca <- prcomp(t(data[rowMeans(data) > 0, ]), scale. = TRUE)
dm  <- destiny::DiffusionMap(pca$x[, 1:10])
dpt <- destiny::DPT(dm, tips = root_index)          # root = a known progenitor cell
# plot DC1/DC2 colored by dpt$dpt / stage / gene (use viridis, not rainbow — sc-conventions)

# Principal-curve lineage fit through an existing embedding (single lineage only).
fit <- princurve::principal_curve(as.matrix(emb[lineage_cells, c("UMAP1","UMAP2")]))
plot(fit); points(emb[, c("UMAP1","UMAP2")])        # fit$lambda = arc-length pseudotime
```

For a more robust R curve/lineage tool use **slingshot** (handles branching). Monocle2/DDRTree is legacy — prefer **monocle3** or slingshot for new work. Diffusion maps need a single true lineage (remove doublets/outliers first) or they distort.
