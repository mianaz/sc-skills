# Figure sizing & typography

House defaults so figures look consistent and publication-ready, even for plot types without a
worked example.

## Typography

| element | value |
|---|---|
| global base | `theme_classic(base_size = 14)` / `theme_bw()` — 14 pt |
| axis text (ticks) | `element_text(colour = "black", size = 12)` |
| axis / plot titles | `size = 16`; center titles `plot.title = element_text(hjust = 0.5)` |
| on-plot cluster labels | `label.size ≈ 4`, bold (`fontface = "bold"`, `family = "sans"`) |
| legend / strip / caption | ~10–12 pt (no fixed rule — keep ≤ axis text) |
| rotated tick labels | 45° (or 60° for long names, 90° for heatmaps), always `hjust = 1` |

Compact multi-panel figures built with cowplot may drop to `theme_cowplot(font_size = 12)` — that's
fine; 14 is the single-panel default, not a hard floor.

## Point size by density (UMAP/embedding)

`pt.size` should shrink as cell count grows. Approximate buckets:

| cells (n) | pt.size |
|---|---|
| dense (≳ 50k) | 0.01 – 0.1 |
| medium (~5–50k) | 0.2 – 0.5 |
| sparse (≲ 5k) | 1 – 1.5 |

Keep vector output (`raster = FALSE`) when feasible; once a vector PDF gets unwieldy (roughly
> 50–100k points) **rasterize the points layer only** (ggrastr, or scop `raster = TRUE`) so axes/text
stay vector. Don't rasterize the whole figure.

## Save dimensions by plot type (inches)

Starting points — scale to the actual panel count / number of categories, don't apply literally.

| plot | width × height |
|---|---|
| single UMAP | 8 × 7 |
| UMAP split into a panel grid | size from panel count (e.g. ~24 × 21 for a large grid) |
| flipped marker dotplot | 10 × ~2.0–2.2 (add height per extra gene-set row) |
| small metric / bar (e.g. LISI) | 5 × 4 |
| trajectory (slingshot/monocle) | 7 × 4.75 |
| heatmap | ~7 × 5 (or cellwidth/cellheight for square tiles) |
| large GO/igraph network | 30 × 30 (`V(g)$label.cex` ~3.8–4.5, bold key nodes) |

Default device is vector **PDF**; also export PNG ≥600 dpi (see scientific-reproducibility → figure-output-contract.md). These sizes are about
keeping the plotting area compact and text legible — pair with the dual-version + source-data rules.
