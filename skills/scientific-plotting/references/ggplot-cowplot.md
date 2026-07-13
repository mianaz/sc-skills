# ggplot2 + cowplot Fallback

Use ggplot2 when tidyplots or scop cannot express the desired figure, or when you need precise control over annotations, facets, scales, and labels. Use cowplot for clean publication styling and plot assembly.

## Base pattern

```r
library(ggplot2)
library(cowplot)

p <- ggplot(df, aes(x, y, color = group)) +
  geom_point(size = 1, alpha = 0.75) +
  scale_color_manual(values = pal) +
  labs(x = "Time (hours)", y = "Response (a.u.)", color = NULL) +
  theme_cowplot(font_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )
```

## Long category labels

```r
p <- ggplot(df, aes(forcats::fct_reorder(pathway, score), score, fill = source)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, 42)) +
  labs(x = NULL, y = "Normalized enrichment score") +
  theme_cowplot(font_size = 11)
```

Use horizontal bars for pathway names, phenotypes, and any long categorical text.

## Repelled labels

```r
p <- ggplot(df, aes(log2_fc, neg_log10_p)) +
  geom_point(aes(color = direction), size = 0.8, alpha = 0.6) +
  ggrepel::geom_text_repel(
    data = subset(df, label_me),
    aes(label = gene),
    box.padding = 0.35,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf,
    size = 3
  ) +
  theme_cowplot(font_size = 12)
```

Prefer fewer intentional labels over many colliding labels.

## P-values with explicit headroom

```r
p <- ggplot(df, aes(group, value, fill = group)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.35) +
  ggbeeswarm::geom_quasirandom(width = 0.12, size = 0.7, alpha = 0.7) +
  ggpubr::stat_compare_means(
    comparisons = list(c("Control", "Treatment")),
    method = "wilcox.test",
    label = "p.format"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none")
```

Always leave enough y-axis expansion for brackets and labels.

## Multi-panel assembly

```r
fig <- cowplot::plot_grid(
  p1, p2, p3,
  labels = c("A", "B", "C"),
  label_fontface = "bold",
  nrow = 1,
  align = "hv",
  axis = "tblr"
)
```

Build panels at compatible sizes before combining. Fix labels and legends at the panel level first.

## Export

```r
readr::write_csv(source_df, "source_data/fig3_source.csv")
ggsave("figures/fig3_clean.pdf", fig, width = 7, height = 3.5, units = "in")
ggsave("figures/fig3_clean.png", fig, width = 7, height = 3.5, units = "in", dpi = 600)
```
