# Ordered-axis trends & gene programs

For data with an **ordered categorical axis** (developmental stages, specification steps, dose) —
not an inferred pseudotime (that's scop `DynamicHeatmap`). Always relevel the axis
to its biological order first. (source: Bian, Gong et al. Nature 2020)

## Single-gene trend: violin + mean-connecting segment

Full distribution at each step plus the directional trend.

```r
library(ggplot2)
ord <- c("CS11","CS12","CS13","CS15","CS17","CS23")            # ordered axis
df  <- data.frame(grp = factor(meta$stage, levels = ord),
                  expr = as.numeric(data[gene, rownames(meta)]))
m   <- tapply(df$expr, df$grp, mean)
ggplot(df, aes(grp, expr, fill = grp)) +
  geom_violin(scale = "width") + geom_point(size = .5) +
  geom_segment(aes(x = seq_along(m)[-length(m)], xend = seq_along(m)[-1],
                   y = m[-length(m)], yend = m[-1]), linewidth = 1, inherit.aes = FALSE) +
  theme_classic() + theme(legend.position = "none") + labs(x = NULL, title = gene)
```
tidyplots alternative for many genes: `add_violin()` + `add_mean_line(group=...)` + `add_sem_ribbon()`,
faceted with `split_plot(by = gene)`.

## Gene-program discovery (ANOVA → PAM → program trends)

Find genes that vary across the ordered axis and group them into co-regulated programs.

```r
# 1. per-stage (or per-cluster) mean matrix, z-scored by gene
avg <- sapply(levels(df$grp), \(g) rowMeans(data[, meta$stage == g, drop = FALSE]))
z   <- t(scale(t(avg)))

# 2. genes varying across the axis: one-way ANOVA, FDR
aov_p <- apply(data[rownames(z), ], 1, \(v)
               summary(aov(v ~ meta$stage))[[1]][["Pr(>F)"]][1])
keep  <- names(which(p.adjust(aov_p, "fdr") < 0.01))

# 3. cluster into k programs (k via factoextra::fviz_nbclust wss + silhouette)
pr <- cluster::pam(z[keep, ], k = 5)$clustering          # or Mfuzz for soft clustering

# 4a. heatmap, genes ordered by program (navy-white-red diverging for z-scores)
pheatmap::pheatmap(z[order(pr), ], cluster_rows = FALSE, cluster_cols = FALSE,
                   annotation_row = data.frame(program = factor(pr[order(pr)])),
                   color = colorRampPalette(c("navy","white","red"))(100), border_color = NA)

# 4b. faceted program trend lines + bootstrapped mean ribbon
long <- reshape2::melt(data.frame(gene = keep, program = factor(pr), z[keep, ], check.names = FALSE),
                       id = c("gene","program"), variable.name = "grp", value.name = "expr")
long$grp <- factor(long$grp, levels = ord)
ggplot(long) +
  geom_line(aes(grp, expr, group = gene, color = program), linewidth = .1) +
  facet_wrap(~program, scales = "free_y") +
  stat_summary(aes(grp, expr, group = 1), fun.data = "mean_cl_boot",
               geom = "smooth", color = "black", alpha = .2, linewidth = 1.1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
```

Each program is a co-regulated module — run GO per program to interpret. Log the clustering choice
(PAM vs Mfuzz, k, z-scoring) in methods.md; it's a judgment call. For genes varying along an inferred
pseudotime (not discrete stages), use scop `RunDynamicFeatures` + `DynamicHeatmap` instead.
