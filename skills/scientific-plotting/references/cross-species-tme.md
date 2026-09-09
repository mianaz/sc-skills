# Cross-species TME figure recipes
(source: sc-paper-distill/papers/2026-humu-tme.md F1–F6 — Courau et al. Nat Immunol 2026)

Use when comparing mouse vs human tumor microenvironments: composition, chemokine source cells,
cell-frequency couplings, NMF gene-program conservation, and coordinated GEP “movements.”
Lock species color: human `#4C78A8`, mouse `#F58518`. Encode species with shape when color is
already used for archetype / genotype.

## Species-split composition violins
(source: papers/2026-humu-tme.md F1)

```r
# (source: sc-paper-distill/papers/2026-humu-tme.md F1)
pal_sp <- c(Human = "#4C78A8", Mouse = "#F58518")
ggplot(freq_df, aes(species, frequency, fill = species)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white") +
  facet_wrap(~ feature, scales = "free_y") +
  scale_fill_manual(values = pal_sp) +
  labs(x = NULL, y = "Frequency") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", strip.background = element_blank())
```

## Cosine similarity of mouse models to human TME archetypes
(source: papers/2026-humu-tme.md F2)

Z-score features on the **pooled** human+mouse matrix, then cosine-similarity. Do not z-score
species separately — that erases the desert vs hot shift.

```r
# (source: sc-paper-distill/papers/2026-humu-tme.md F2)
library(ComplexHeatmap); library(circlize)
X <- scale(feature_mat)                 # samples x features
Xn <- X / sqrt(rowSums(X^2))
S <- tcrossprod(Xn)
Heatmap(S[mouse_ids, human_ids], name = "cosine",
        col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
        rect_gp = gpar(col = "white", lwd = 0.4),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2")
```

On a composition UMAP, human = circles, mouse = triangles.

## Compartment-scaled chemokine heatmap
(source: papers/2026-humu-tme.md F3)

Z-score each gene **across compartments within species**, then column-split human | mouse.
This shows which cell type owns the chemokine, not absolute expression.

```r
# (source: sc-paper-distill/papers/2026-humu-tme.md F3)
scale_within <- function(mat) t(scale(t(mat)))
Heatmap(cbind(scale_within(human_cpmt), scale_within(mouse_cpmt)),
        name = "z", col = colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426")),
        column_split = factor(c(rep("Human", ncol(human_cpmt)),
                                rep("Mouse", ncol(mouse_cpmt)))),
        cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = 0.3))
```

## Human-ordered dual correlation matrices
(source: papers/2026-humu-tme.md F4)

Cluster the human Pearson matrix; plot the mouse matrix in that same order. Discordant pairs
get a sensitivity scatter restricted to desert / mouse-like human samples.

```r
# (source: sc-paper-distill/papers/2026-humu-tme.md F4)
Rh <- cor(human_freq, use = "pairwise.complete.obs")
ord <- hclust(as.dist(1 - Rh))$order
Rm <- cor(mouse_freq, use = "pairwise.complete.obs")
corrplot::corrplot(Rh[ord, ord], method = "color", tl.col = "black", tl.cex = 0.7)
corrplot::corrplot(Rm[ord, ord], method = "color", tl.col = "black", tl.cex = 0.7)
```

## NMF gene-weight scatter
(source: papers/2026-humu-tme.md F5)

Jaccard of top 20 (T) / top 50 (myeloid) genes is the match; the loading scatter is the QC.
Do not call a program conserved from Jaccard > 0.05 alone.

```r
# (source: sc-paper-distill/papers/2026-humu-tme.md F5)
cutoff <- 40
ggplot(w, aes(human_weight, mouse_weight)) +
  geom_point(aes(color = rank_min <= cutoff), size = 1.2, alpha = 0.85) +
  geom_hline(yintercept = sort(w$mouse_weight, decreasing = TRUE)[cutoff],
             linetype = 2, color = "grey50") +
  geom_vline(xintercept = sort(w$human_weight, decreasing = TRUE)[cutoff],
             linetype = 2, color = "grey50") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey70")) +
  coord_equal() + theme_classic(base_size = 11)
```

## 2×2 Kaplan–Meier for a T × myeloid GEP movement
(source: papers/2026-humu-tme.md F6)

Score T-cell cytotoxicity and myeloid IFN programs separately, split each at the cohort median,
plot four curves. Do not collapse to a single myeloid-high vs low split.

```r
# (source: sc-paper-distill/papers/2026-humu-tme.md F6)
library(survival); library(survminer)
clin$T3  <- ifelse(clin$T3_score  >= median(clin$T3_score),  "Hi", "Lo")
clin$My2 <- ifelse(clin$My2_score >= median(clin$My2_score), "Hi", "Lo")
clin$axis <- factor(paste0("T3", clin$T3, "_My2", clin$My2),
                    levels = c("T3Hi_My2Hi", "T3Hi_My2Lo", "T3Lo_My2Hi", "T3Lo_My2Lo"))
ggsurvplot(survfit(Surv(time, event) ~ axis, data = clin), pval = TRUE, risk.table = TRUE,
           palette = c("#1B9E77", "#7570B3", "#D95F02", "#666666"))
```
