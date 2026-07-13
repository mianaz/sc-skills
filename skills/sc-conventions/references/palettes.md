# Palettes & color conventions

Source this once per session so every plot is consistent.

## Categorical (discrete groups)
```r
library(RColorBrewer)
# "Pair2"-style: the 12-color Paired set
pal_categorical <- function(n) {
  base <- RColorBrewer::brewer.pal(12, "Paired")
  grDevices::colorRampPalette(base)(n)
}
# ggpubr palettes are also acceptable, e.g. ggpubr::get_palette("npg", n)
```

## Many categories (high cardinality, >12 groups)

The Paired ramp above interpolates new hues between only 12 anchors, so at ~20+ groups
(e.g. a 38-cluster UMAP) adjacent clusters become near-indistinguishable. For high
cardinality use a palette engineered for maximal categorical separability instead.

```r
# Option A — perceptually-spaced (PREFERRED for publication): Polychrome / pals.
pal_many <- function(n) {
  if (requireNamespace("Polychrome", quietly = TRUE)) {
    seed <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B","#E377C2","#17BECF")
    return(unname(Polychrome::createPalette(n, seedcolors = seed, range = c(20, 90))))
  }
  if (requireNamespace("pals", quietly = TRUE))
    return(unname(if (n <= 36) pals::polychrome(n) else pals::glasbey(n)))
  pal_categorical(n)                       # fall back to the Paired ramp
}
# NOTE: Polychrome::createPalette draws stochastically — set a seed for reproducibility.

# Option B — concatenate every Brewer qualitative set, then recycle (no extra pkgs).
# (source: Xue et al. Nature 2022)
pal_many_brewer <- function(n) {
  qual <- subset(RColorBrewer::brewer.pal.info, category == "qual")
  cols <- unlist(lapply(rownames(qual), \(p) RColorBrewer::brewer.pal(qual[p, "maxcolors"], p)))
  rep(cols, length.out = n)
}
```

Bind colors to names (`names(cols) <- levels(factor(seu$cluster))`) so the mapping stays
consistent across every panel of a figure. On dense high-k UMAPs, drop on-plot labels
(use a multi-column legend) and shrink `pt.size`. Prefer Option A for manuscripts;
Option B is the dependency-free fallback (its hues are distinct but not perceptually tuned).

## Same color AND same ORDER across panels (UMAP legend ↔ dotplot)

A category must appear in the **same order** in every panel, not just the same color.
The classic bug: the marker dotplot is ordered one way (e.g. alphabetical) while the
UMAP legend is ordered another (e.g. by cluster frequency), so a reader cannot line up
"the green one" between panels. Define **one canonical order** (a `levels` vector) and
apply it to BOTH the dotplot y-axis/panels and the UMAP factor.

```r
ct_lv <- canonical_order[canonical_order %in% present_labels]  # one order, everywhere
seu$celltype <- factor(seu$celltype, levels = ct_lv)           # UMAP legend order
# dotplot: order y + marker panels by the same ct_lv
```

### Gotcha: `scop::CellDimPlot` applies `palcolor` POSITIONALLY, ignoring names

`scop::CellDimPlot` orders groups by **frequency** and maps `palcolor` by **position**,
so a *named* palette vector does NOT bind by name — e.g. the largest cluster silently
takes `palcolor[1]` regardless of its label. On a plain **character** column this gives
the UMAP a different color/order per type than a name-bound dotplot (Fibroblast drawn
teal instead of red). Fix: make the group column a **factor with explicit levels** and
pass `palcolor` **in that same level order**:

```r
seu$celltype <- factor(seu$celltype, levels = ct_lv)
scop::CellDimPlot(seu, group.by = "celltype",
                  palcolor = my_named_pal[ct_lv])   # vector in factor-level order
```
Then the UMAP legend order + colors match the name-bound dotplot exactly. (Plain ggplot
`scale_*_manual(values = named_vec)` DOES respect names — the positional quirk is
scop-specific; a login-safe ggplot UMAP from the embedding CSV is a good cross-check.)

### Gotcha: a cowplot/patchwork-composed dotplot title can't be stripped afterward

If a grouped dotplot bakes its title into a sub-panel (e.g. a `cowplot::plot_grid`
bracket strip), `+ ggtitle(NULL)` on the returned object adds a blank title to the
OUTER canvas and leaves the real one — so the "clean" version is byte-identical to the
"titled" one. Render the clean version from scratch with `title = NULL` instead of
trying to remove it.

## Continuous / diverging (heatmaps, signed scores) — decoupleR style
```r
colors      <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use  <- grDevices::colorRampPalette(colors = colors)(100)
my_breaks   <- c(seq(-1.25, 0, length.out = ceiling(100 / 2) + 1),
                 seq(0.05, 1.25, length.out = floor(100 / 2)))
```

## Sequential / continuous (non-negative: expression, pseudotime, density)

For 0→max data (gene expression on feature plots, pseudotime, density) use a **perceptually-uniform
sequential** palette — viridis family. Reserve the RdBu diverging ramp above for *signed* data centered
at 0. (source: Bian, Gong et al. Nature 2020)

```r
viridisLite::viridis(100)                      # or "magma" / "inferno" / "plasma" / "cividis"
# two-tone for marker FeaturePlots (zero = light grey): c("#E6E6E6", "red")
# clamp z-scores before a heatmap so a few extreme cells don't wash out the scale:
z[z >  1.5] <-  1.5; z[z < -1.5] <- -1.5
```

**Avoid rainbow/jet** (`blue→green→yellow→red`): perceptually non-uniform and not colorblind-safe,
even though it is common in published figures. Default to viridis. This satisfies the scientific-plotting
colorblind-friendly non-negotiable.

## Heatmap visual style (square tiles, white borders)
```r
pheatmap::pheatmap(
  mat            = mat,
  color          = colors.use,
  breaks         = my_breaks,
  border_color   = "white",
  cellwidth      = 20,
  cellheight     = 20,
  treeheight_row = 20,
  treeheight_col = 20
)
```
Adjust the symmetric `my_breaks` range to the data's score range; keep it centered at 0.
