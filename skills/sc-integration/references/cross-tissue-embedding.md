# Cross-tissue embedding of a broadly distributed compartment

To study **one cell class across many tissues at once** (blood / endothelial / epithelial / myeloid
across every organ), don't just merge everything and hope — the largest tissue dominates the embedding
and swamps the shared program you're trying to see. Downsample per type-per-tissue, then build the
embedding on a **marker-gene set**, not all HVGs.
(source: https://doi.org/10.1126/science.aba7721)

## Recipe

1. **Extract** the annotated cells of the class from every tissue.
2. **Downsample to ≤N per (cell type × tissue)** — Cao used ≤5000; take all cells for types with fewer.
   *Why:* equal footing per group so a 2M-cell organ can't define the axes; the rare-tissue signal survives.
3. **Feature set = top per-type markers**, not global HVGs: per type, keep genes with `q < 0.05` and
   ≥2× fold vs the second-ranked type, ordered by median q across tissues.
   *Why:* HVGs re-encode the tissue-vs-tissue and technical variance you want to look *past*; a
   type-marker set focuses the embedding on cell-identity structure.
4. **Embed:** PCA on that gene set → UMAP (`n_neighbors≈50, min_dist=0.1, metric="cosine"`) →
   Louvain at low resolution (Cao: `res≈1e-4`) → annotate clusters by tissue-of-origin + markers.

```r
# Seurat sketch (Monocle3 in the paper; the choices transfer)
sub <- subset(seu, subset = compartment == "endothelial")
# per (cell_type, tissue) cap
cells <- unlist(lapply(split(colnames(sub), paste(sub$cell_type, sub$tissue)),
                       function(cc) if (length(cc) > 5000) sample(cc, 5000) else cc))
sub <- sub[, cells]
marker_genes <- top_type_markers(sub, group.by = "cell_type", q = 0.05, min_fc = 2)   # your DE helper
sub <- Seurat::ScaleData(sub, features = marker_genes) |>
       Seurat::RunPCA(features = marker_genes) |>
       Seurat::RunUMAP(dims = 1:30, n.neighbors = 50, min.dist = 0.1, metric = "cosine")
```

- **When to reach for it:** cross-organ / cross-region atlases where you want organ-specific
  *specializations* of a shared type (e.g. tissue-resident macrophage flavours), not a per-organ view.
- To bring in an **external reference atlas** of the same compartment, co-embed with Seurat anchors on
  shared HVGs (see harmony.md / scvi-scanvi.md for the batch-correction layer).
- Log the per-group cap, the marker-set definition, and the resolution in methods.md.
