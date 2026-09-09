---
name: scientific-plotting
description: Use when making, revising, or debugging an R scientific figure (ggplot, scop, heatmap, volcano, composition). Prefer tidyplots, scop, or ggplot2+cowplot.
---
# Scientific Plotting

## Scope

Use this as the R-first plotting skill for scientific figures, from quick analysis plots to manuscript panels. It is a general scientific plotting skill (formerly `sc-plotting`), not only single-cell.

For Python-only matplotlib/seaborn work or journal-template sizing, route to `scientific-visualization`. For raw statistical modeling, use the relevant stats skill first, then return here for the figure. For designing a color *system* for a web/artifact/dashboard context, that is `dataviz`; here palettes are fixed (Paired categorical, decoupleR RdBu continuous) per `sc-conventions` and kept semantically consistent across panels — do not re-derive them.

## Workflow

1. **Contract first (manuscript/multi-panel figures only).** Write the one-sentence claim the figure must defend, then map each planned panel to that claim and drop any panel that carries no unique piece of evidence. Skip this for standalone exploratory/analysis plots — go straight to step 2. (adapted from nature-figure)
2. Classify the data and figure: tidy dataframe, Seurat/scop object, matrix/heatmap, DE/enrichment result, spatial object, time course, or custom layout.
3. Read the matching reference before coding. Do not rely on memory for package-specific arguments.
4. Build from a tidy source table or an explicit matrix. Save the table/matrix as figure source data.
5. Plot the clean version first, then add the titled/annotated version. Keep both.
6. Run a legibility pass: labels cannot overlap, axis text must fit, p-values must be readable, and legends must not dominate the data.
7. Export vector PDF plus high-resolution PNG/TIFF where useful. Use at least 600 dpi for line-art PNGs. Keep figure text editable (`references/stats-and-integrity.md`).

## Tool Router

| Need | First choice | Reference |
|---|---|---|
| Tidy dataframe plots: boxes, violins, beeswarms, bars, time courses, stacks, heatmaps, exact p-values | tidyplots | `references/tidyplots.md` |
| Label overlap, long categories, cramped axes, p-value brackets, legend placement | layout checklist | `references/layout-legibility.md` |
| Seurat/single-cell/spatial core figures: UMAP, feature maps, composition, marker dotplots, density, group heatmaps | scop | `references/scop-plots.md` |
| scop specialized plots: DE, enrichment/GSEA/GSVA, cell-cell communication, deconvolution, spatial gradients, metabolism, SCENIC, ATAC tracks | scop | `references/scop-specialized.md` |
| Custom ggplot, annotations, facets, patchwork/cowplot, label repulsion | ggplot2 + cowplot | `references/ggplot-cowplot.md` |
| Colored facet strips, marker/GSEA dotplots with grouped strips | ggh4x | `references/ggh4x-facets.md` |
| Gene trends across ordered stages/steps + gene-program (module) discovery & trend panels | ggplot/pheatmap | `references/trends-and-programs.md` |
| Volcano plots from DE tables | scop or tidyplots/ggplot fallback | `references/volcano.md` |
| Simple matrix heatmap | pheatmap/RdBu | `references/heatmaps.md` |
| Annotated heatmap with row/column annotations, splits, side bars | ComplexHeatmap | `references/complexheatmap.md` |
| Human vs mouse TME: composition violins, archetype cosine heatmaps, chemokine source-cell heatmaps, frequency-coupling matrices, NMF gene-weight scatters, 2×2 GEP-movement KM | ggplot2 + ComplexHeatmap | `references/cross-species-tme.md` |
| Per-panel stats reporting, editable-text export, image-integrity checklist | reporting/QA contract | `references/stats-and-integrity.md` |

## Non-Negotiables

- Never leave overlapping labels, clipped labels, unreadable rotated text, or p-value brackets running into points.
- Prefer exact p-values over stars unless the target venue or user explicitly wants stars.
- Use colorblind-friendly palettes and preserve semantic colors across related panels.
- Show individual observations for small to medium n; summarize with mean/median/error only as an overlay.
- For many categories, increase figure width, wrap labels, facet, flip coordinates, or move labels to a source-data table rather than shrinking text into illegibility.
- Document the statistical test and adjustment method in code comments, methods text, or the figure caption.
- **Record stats and provenance per panel** with the structured block in `references/stats-and-integrity.md` (n, replicate definitions, center, spread, test, correction, p display, source-data file); for image panels record the image-integrity block there too.
- **Data fidelity:** a title/subtitle that states a claim ("X higher in treated") must be verifiable against every plotted row; never let excluded/filtered data silently enter a summary, and use one canonical number for a given quantity across all panels and captions.
- **Render-then-verify:** after saving, look at the actual output (open the PNG / crop panels) and check for label/tick collisions, leader-line crossings, and colour confusion before calling a figure done — do not trust the code alone. (Both adapted from Anthropic Claude Science `figure-style`.)
