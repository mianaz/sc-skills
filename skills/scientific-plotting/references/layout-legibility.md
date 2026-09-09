# Layout and Legibility Checklist

Use this before considering a plot finished. The most common failure is not the plot type; it is unreadable labels.

## Axis labels

- Keep horizontal labels only when they are short and fewer than about 6-8 categories.
- Rotate x-axis labels to 45 degrees for moderate labels; set `hjust = 1`.
- Use 90 degrees mainly for heatmaps or very short categorical labels. If 90-degree labels are still cramped, widen the figure or wrap.
- Wrap long labels before plotting: `stringr::str_wrap(label, width = 18)` for compact axes, `width = 25-40` for pathway names.
- Prefer `coord_flip()` for ranked bars, pathway names, or any category labels longer than a few words.
- Use explicit ordering with `forcats::fct_reorder()`, `factor(levels = ...)`, or tidyplots `reorder_x_axis_levels()`.

## Label overlap

Use repulsion or reduce the number of labels; never keep colliding labels.

```r
library(ggrepel)
p + ggrepel::geom_text_repel(
  aes(label = label),
  box.padding = 0.35,
  point.padding = 0.2,
  max.overlaps = Inf,
  min.segment.length = 0
)
```

For tidyplots:

```r
p |>
  add_data_labels_repel(
    label = gene,
    data = filter_rows(label_me),
    background = TRUE,
    box.padding = 0.25,
    max.overlaps = Inf
  )
```

For dense UMAPs, label broad cell types only. If there are more than 12-15 groups, prefer a legend or split panels over `label = TRUE`.

### Indexed labels for high-cardinality embeddings

When even a legend with 20-40 long names is hard to map back to the plot, **number each
population at its median** and decode the numbers with a separate key. Clean on the plot,
and the key can be a standalone vector PDF for panel assembly (see scientific-reproducibility → figure-output-contract.md).
(source: sc-paper-distill/papers/2019-FCA-liver.md F1, F2 — Popescu et al. Nature 2019)

```r
library(ggplot2); library(dplyr)
df  <- data.frame(DR1 = emb[, 1], DR2 = emb[, 2], lab = factor(labels))
med <- df |> group_by(lab) |> summarise(DR1 = median(DR1), DR2 = median(DR2)) |>
  mutate(idx = row_number())
ggplot(df, aes(DR1, DR2, color = lab)) +
  geom_point(size = 0.2) +
  geom_point(data = med, shape = 21, fill = "gray", color = "gray30", size = 6, alpha = .5) +
  geom_text(data = med, aes(label = idx), color = "black", size = 3) +
  scale_color_manual(values = pal_many(nlevels(df$lab))) +   # sc-conventions pal_many
  theme_classic() + theme(legend.position = "none")          # key rendered separately
```

scop `CellDimPlot` remains the default for ≤~15 groups; this is the numbered-index variant
for very high cardinality.

### Highlight a subset on an embedding (grey the rest)

To emphasize one stage/group/lineage on a UMAP, grey out all other cells **and draw the highlighted
cells last** so they aren't hidden under the grey context layer. Pair with a composition pie if the
quantitative makeup matters. (source: sc-paper-distill/papers/2020-human-macrophage-dev.md F3 — Bian, Gong et al. Nature 2020)

```r
df$grp <- ifelse(df$stage == "CS15", as.character(df$cluster), "others")
df <- df[order(df$grp == "others", decreasing = TRUE), ]   # 'others' first → drawn underneath
ggplot(df, aes(UMAP1, UMAP2, color = grp)) + geom_point(size = .8) +
  scale_color_manual(values = c(pal_many(n), others = "grey85")) + theme_bw()
# composition pie: ggpubr::ggpie(counts, "n", fill="cluster") + ggrepel::geom_label_repel(...)
```
Use a light grey (`grey85`/`#E6E6E6`) for context, consistent with the volcano background convention.

## P-values and annotations

- Use exact p-values by default: `p = 0.013`, not only `*`.
- If there are multiple comparisons, adjust p-values and state the method.
- If brackets overlap data, add y-axis padding or move the test text outside the data cloud.
- In tidyplots, use `padding_top`, `step.increase`, and `bracket.nudge.y` in `add_test_pvalue()`.
- In ggplot2, reserve vertical space with `scale_y_continuous(expand = expansion(mult = c(0.05, 0.18)))`.

## Legends

- Put legends on the right for tall plots and on the bottom for wide plots.
- For many groups, use multiple columns: `guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))`.
- Remove redundant legends when direct labels or facet labels already identify groups.
- Do not let a legend consume more visual weight than the plotted data.

## Figure size

- Single panel with short labels: 3.5 x 3 in is often enough.
- Discrete x-axis with 8-15 labels: start around 5-7 x 3.5 in.
- Long pathway/category labels: use horizontal bars, 5-7 x 4-6 in.
- Multi-panel/faceted plots: set dimensions from the number of panels, not from a default device.
- For tidyplots, use `adjust_size(width = ..., height = ..., unit = "mm")`.
- For ggplot2, use `ggsave(width = ..., height = ..., units = "in", dpi = 600)`.

## Multi-panel composition

Individual, standalone figures are first-class. They are what you reuse in slides,
posters, and for quick customization, and most analysis outputs should be a single
clean figure. Compose a multi-panel figure only when the panels form one argument (a
manuscript display item). Both are valid end products — never force a standalone plot
into a grid just to fill a canvas.

When you *do* assemble a panel:

- Classify it as one of: quantitative grid, schematic-led composite, image plate + quant,
  or asymmetric mixed-modality.
- Prefer **one hero panel plus subordinate evidence panels** over a uniform grid of
  equal-sized subplots. Size panels by importance, not by symmetry.
- Give each panel a unique job; progress overview → deviation → relationship so no two
  panels carry the same evidence.
- Keep one restrained palette across the whole panel and reuse it semantically across
  sub-panels (sc-conventions palettes); shared or direct labels beat one legend per
  sub-panel.
- Panel tags: lowercase bold letters near the top-left, not large badges.

(adapted from nature-figure archetypes / 2026 Nature sample observations)

## Final Audit

Before saving final outputs, check:

- Axis titles include units where applicable.
- Tick labels fit and are not clipped.
- Every color mapping has a readable legend or direct label.
- Text is readable at final intended size.
- No panel title, strip label, p-value, or legend overlaps data.
- Source data for the figure has been saved.
