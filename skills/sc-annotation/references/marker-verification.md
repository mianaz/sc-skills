# Canonical-marker verification (MANDATORY gate)

No label is trusted until confirmed by canonical markers with citations.

1. Pull candidate cell types' markers + PMID/DOI from references/markers.md.
2. Build a marker **dotplot** via `scop::GroupHeatmap(features, group.by, add_dot = TRUE)`
   (genes × cell types; see memory `scop-marker-dotplot.md`). Do **not** use
   `FeatureStatPlot(plot_type = "dot")` (facets one gene per panel). Genes grouped by
   cell type on one axis, clusters/types on the other.
3. Confirm each cluster expresses its assigned type's canonical markers.
4. On disagreement between automated calls and marker evidence, favor the
   (documented) marker evidence.

## Two-level workflow (default)
1. Coarse compartments (T/NK, myeloid, B/plasma, stromal, endothelial, epithelial) → verify.
2. Subcluster within each compartment.
3. Fine subtypes → verify again.

Record per-cluster: assigned identity, supporting markers + citations, and whether
automated calls agreed, in methods.md.
