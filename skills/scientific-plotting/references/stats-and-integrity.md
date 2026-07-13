# Stats Reporting, Editable-Text Export, and Image Integrity

The reporting and export contract a publication-ready figure must satisfy. Load when
finalizing a manuscript/multi-panel figure, before delivery, or whenever a figure
contains microscopy/spatial images, blots/gels, or a statistical claim.
(adapted from nature-figure / Yuan1z0825/nature-skills; always verify the target
journal's current author guide before submission.)

## Per-panel statistics block

For every quantitative panel, record — in the caption, methods, or a sidecar file:

```text
n definition:
biological replicates:
technical replicates:
center statistic:                 # mean / median
spread / interval:                # SD / SEM / 95% CI / IQR
test:
multiple-comparison correction:
p-value display:                  # exact p preferred over stars
source-data file:
```

For model / machine-learning panels, also record:

```text
train/validation/test split:
number of seeds or folds:
metric definition:
confidence-interval or variability definition:
baseline definition:
```

This is the structured form of the "document the statistical test" non-negotiable — one
filled block per panel, not a single line for the whole figure.

## Editable-text, font-embedded export (R)

Vector output must keep text **selectable** so a copyeditor can restyle labels without
re-rendering. Size in mm (journals think in mm); convert with `mm / 25.4`.

```r
w_mm <- 89; h_mm <- 70            # 89 mm ≈ single-column, 183 mm ≈ double-column

svglite::svglite("fig.svg", width = w_mm/25.4, height = h_mm/25.4)      # editable text
print(p); dev.off()

grDevices::cairo_pdf("fig.pdf", width = w_mm/25.4, height = h_mm/25.4, family = "Arial")
print(p); dev.off()

ragg::agg_tiff("fig.tiff", width = w_mm/25.4, height = h_mm/25.4, units = "in", res = 600)
print(p); dev.off()
```

Then **open the exported SVG/PDF** and confirm: text is selectable (not outlined),
labels do not overlap, and the figure still reads at final printed size. This is the
export half of render-then-verify — trusting the plotting code is not enough.

For journal-template sizing beyond this (per-journal column widths, height caps, full
spec tables), route to `scientific-visualization`; do not duplicate those tables here.

## Image-integrity block (microscopy / spatial / blots)

For each image panel, record:

```text
raw file:
processed file:
crop:
brightness / contrast / gamma:
pseudo-color:
scale-bar calibration:            # calibrated units, not just a magnification factor
stitching:
reuse in other figures:
quantification link:
```

Global adjustments (whole-image brightness/contrast) are safer than local selective
edits. If an adjustment changes the visibility of relevant background or signal, flag it
rather than silently normalizing it away. Every image panel needs a calibrated scale bar.
