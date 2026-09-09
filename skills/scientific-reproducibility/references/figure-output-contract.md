# Figure output contract (apply to every reported figure)

Language-agnostic in intent; R helpers below, Python equivalents noted.

## Dual-version (always save both)
- `_titled`: title + subtitle describing the plot and the full stat-test details (test name, comparison, n, exact p-value).
- `_clean`: only necessary components; NO title/subtitle unless the figure is unintelligible without them.

## Source data (always save)
Write the tidy data behind the plot to `source_data/<plot>.csv` (one csv per plot) and the plot object to `source_data/<plot>.rds`. This lets the figure be re-styled later or rebuilt elsewhere (e.g. GraphPad Prism). In Python, save the tidy `DataFrame` to csv and pickle the figure/axes or the plotting inputs.

## Statistics
Wherever a comparison is shown, run the appropriate test and display the **exact p-value**, not significance stars (unless stars are explicitly requested). Pick the test by data type (e.g. Wilcoxon for two non-normal groups, Kruskal–Wallis + Dunn for >2, proportion tests for compositions) and record the choice in `methods.md`.

## State exclusions on the figure
When a panel drops observations before plotting — especially **composition / proportion / box plots**
where a low-n group would read as a real zero — state the exclusion rule in the caption or subtitle
(e.g. "samples with ≤200 blood cells excluded") and keep the dropped rows flagged in the source-data csv.
(source: sc-paper-distill/papers/2020-descartes-fetal-atlas.md §4 — Cao et al. Science 2020)

- **Why:** an unstated threshold makes a filtered plot indistinguishable from complete data — a reader
  can't tell an excluded group from a genuinely absent one. This is the concrete form of the data-fidelity
  rule (never let excluded data silently enter a summary): the exclusion must be *visible*, not just applied.
- Record the same threshold + count dropped in `methods.md`; the csv keeps every row with an `included`
  column so the figure is rebuildable and the exclusion auditable.

## Sizing & resolution
- Constrain the plotting area so the plot is compact (not visually "empty").
- Labels, axis titles, and tick labels large and clearly legible.
- Save PNG at ≥300 dpi (prefer 600) AND a vector PDF. (Python: `savefig(..., dpi=600)` + a `.pdf`.)

## Separate legend for multi-panel assembly
When a figure is assembled from several panels (Illustrator/cowplot/patchwork), export the legend as its **own vector PDF** and the panels legend-stripped, so a shared key appears once and panels rearrange freely.

```r
p_clean    <- p + ggplot2::theme(legend.position = "none")          # panel, no legend
legend_pdf <- cowplot::get_legend(p + ggplot2::theme(legend.box.margin = ggplot2::margin(0,0,0,6)))
ggplot2::ggsave("figures/<plot>_clean.pdf", p_clean, width = w, height = h)
ggplot2::ggsave("figures/<plot>_legend.pdf", cowplot::ggdraw(legend_pdf), width = 2, height = 3)
```
For a numbered/indexed key (high-cardinality categories), render it as a standalone grid of numbered swatches (`theme_void`). For ComplexHeatmap, `draw(ht, show_heatmap_legend = FALSE)` for the body + a second PDF with just the `Legend()`.

## Canonical save helpers
```r
save_plot_pair <- function(p_titled, p_clean, basename, width, height,
                           dir_fig = "figures") {
  for (v in list(list("titled", p_titled), list("clean", p_clean))) {
    tag <- v[[1]]; p <- v[[2]]
    ggplot2::ggsave(file.path(dir_fig, sprintf("%s_%s.png", basename, tag)),
                    p, width = width, height = height, dpi = 600)
    ggplot2::ggsave(file.path(dir_fig, sprintf("%s_%s.pdf", basename, tag)),
                    p, width = width, height = height)
  }
}
save_source_data <- function(tidy_df, plot_obj, basename,
                             dir_src = "source_data") {
  utils::write.csv(tidy_df, file.path(dir_src, paste0(basename, ".csv")),
                   row.names = FALSE)
  saveRDS(plot_obj, file.path(dir_src, paste0(basename, ".rds")))
}
```
(For non-ggplot outputs like pheatmap, use the appropriate device: `png(..., res = 600)` / `pdf(...)`.)

Use a single `<plot>` basename for all four artifacts (titled/clean × png/pdf) plus the csv/rds, so a figure and its source data are trivially linked.
