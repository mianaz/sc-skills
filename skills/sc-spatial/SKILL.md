---
name: sc-spatial
description: Use when analysing spatially resolved transcriptomics (Visium, MERFISH, Slide-seq, Xenium, MERSCOPE).
---

# Spatial Transcriptomics

## When to use
- User has Visium, MERFISH, Slide-seq, Xenium, or similar spatially resolved data
- User wants spatially variable genes, tissue domains/niches, or neighbourhood enrichment
- User wants to overlay gene expression on tissue H&E images

## When NOT to use
- Standard scRNA-seq with no spatial coordinates → `sc-preprocessing` + `sc-annotation`
- Spot resolution too coarse to resolve individual cells without deconvolution → RCTD or cell2location first

## Method decision tree
1. Load data: `sq.read.visium()` (Visium), or `anndata.read_h5ad()` for pre-processed objects.
2. QC and normalisation: scanpy QC (filter low-count spots) → `sc.pp.normalize_total` → `sc.pp.log1p`.
3. Spatially variable genes: `sq.gr.spatial_autocorr(mode="moran")` — higher I = more spatially structured expression.
4. Spatial clustering: `sq.gr.spatial_neighbors` → Leiden on expression + spatial adjacency.
5. Neighbourhood enrichment: `sq.gr.nhood_enrichment` — needs cell-type labels (deconvolution or sub-spot platforms like Xenium).
6. Visualisation: `sq.pl.spatial_scatter` or `sc.pl.spatial`.

## Key caveats
- **Spot ≠ cell (Visium):** each 55 µm spot captures ~2–10 cells; run RCTD or cell2location for per-spot cell-type proportions.
- Moran's I is spatial autocorrelation, not causation.
- Xenium/MERFISH: cell segmentation quality is the dominant QC concern.
- Do not pool tissue sections without batch correction.

## Recipes (load when needed)
- Two-type Visium hex-grid colocalization → [references/visium-colocalization.md](references/visium-colocalization.md)
- H&E segmentation → spot-to-region assignment → [references/he-region-assignment.md](references/he-region-assignment.md)

## Routes to
- `scanpy` — normalisation, clustering, AnnData
