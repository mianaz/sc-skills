---
paper: "human-macrophage-dev — Deciphering human macrophage development at single-cell resolution"
authors: "Bian Z, Gong Y, et al."
doi: "10.1038/s41586-020-2316-7"
journal: "Nature"
year: 2020
repo: "https://github.com/yandgong307/human_macrophage_project @ 5efb6a623d8f69a17872d44a2bfbd1d441e61220"
data_accession: "GEO/GSA — see paper (STRT-seq + 10x scRNA-seq of human embryonic/fetal tissues)"
distilled: 2026-06-25
tags: [human, embryo, yolk-sac, fetal-liver, macrophage, microglia, development, haematopoiesis, EMP, STRT-seq, 10x, trajectory, SCENIC]
---

# human-macrophage-dev — distillation

## 1. Overview
- **What it is & why notable:** Maps human macrophage development from yolk-sac EMPs/HSPCs across embryonic tissues (yolk sac, liver, head/microglia, skin, lung, blood) using plate-based **STRT-seq** + droplet **10x**. Chosen for **method breadth on developmental/time-course data**: ordered-stage gene-pattern discovery, several R-native trajectory methods (monocle2, destiny diffusion + DPT, principal curve), cross-platform integration, and a set of developmental-biology plot idioms (violin-with-mean-trend, stage-highlight UMAP + pie, circular dendrogram).
- **Data:** human; embryonic/fetal yolk sac, liver, head, skin, lung, blood; Carnegie stages CS11–CS23; STRT-seq (full-length, plate) + 10x; integrated with external adult datasets for comparison.
- **Code availability:** ships analysis + figure code (R/Seurat **v2** API — `@dr`, `FindVariableGenes`, `SetAllIdent`, `do.identify`; modernize to v5 on reuse). Read at snippet fidelity: `STRT-seq_analysis.r` (1211 lines — main analysis + Fig1/3/4/7 blocks), `data integration of STRT-seq and 10x.R`, `data integration of cs11 and cs15 in YS.R`. **Missing / external** (sourced but not in repo): `plot_gene_violin_matrix.R`, `get_smoothdata.R`, `go_analysis.R`, `plot_mpath.R`, and the SCENIC run (only its precomputed `binaryRegulonActivity` Rds is consumed). Paper full text is paywalled (no OA); figure mapping uses the script's own `Fig` section comments.

## 2. Methods catalog

### M1. Ordered-stage gene-pattern discovery (ANOVA → PAM clustering → programs)
- **Purpose:** find genes that vary across an ordered axis (developmental stages or specification clusters) and group them into co-regulated temporal programs.
- **How used:** per-gene one-way **ANOVA** across groups (FDR-adjusted) on cluster/stage-averaged expression → z-score the averaged matrix → choose k by `fviz_nbclust` (wss + silhouette) → **PAM** (`cluster::pam`) into k programs → order genes by program → pheatmap + faceted trend lines + GO per program.
- **Tools / packages:** stats `aov`, `cluster::pam`, `factoextra::fviz_nbclust`, pheatmap.
- **Code:** `STRT-seq_analysis.r:691-757` (anova.test + pam + line plots), `:834-848` (top TF heatmap), `:949-999` (liver version).
- **Interpretation:** each program = a temporal/lineage gene module (e.g. genes up across head-macrophage maturation), interpretable via GO.
- **Reusable?** yes → scientific-plotting (program heatmap + faceted trend lines) + the anova→pam recipe.

### M2. R-native diffusion pseudotime (destiny DiffusionMap + DPT)
- **Purpose:** order a single macrophage lineage and visualize maturation in diffusion space, in R.
- **How used:** `prcomp` (scaled) → `destiny::DiffusionMap(pca$x[,1:10..15])` → `DPT(dm, tips=)`; plot DC1/DC2 colored by stage / DPT / gene; inspect `eigenvalues(dm)`.
- **Tools / packages:** destiny, ggbeeswarm (`geom_quasirandom` of DC1/DPT by stage).
- **Code:** `STRT-seq_analysis.r:872-947` (liver), `:1001-1048` (skin).
- **Interpretation:** R-native equivalent of scanpy DPT; the suite currently routes diffusion pseudotime to scanpy only.
- **Reusable?** yes → sc-trajectory (R/destiny diffusion-pseudotime path).

### M3. Principal-curve lineage fit
- **Purpose:** fit a smooth 1D lineage curve through cells in an existing embedding (lightweight trajectory).
- **How used:** `princurve::principal_curve(embedding[lineage_cells, ])`; overlay with base `plot(fit)` + `whiskers()` showing each cell's projection.
- **Code:** `STRT-seq_analysis.r:499-522`.
- **Reusable?** yes → sc-trajectory (princurve as a simple lineage-fit option alongside PAGA/monocle/DPT).

### M4. Monocle2 DDRTree trajectory (wrapped)
- **Purpose:** branched pseudotime of myeloid specification.
- **How used:** `do_monocle()` wrapper — `newCellDataSet(negbinomial)` → dispersion-based ordering genes (minus cell-cycle) → `reduceDimension(DDRTree, param.gamma)` → `orderCells` (optional root) → `plot_cell_trajectory`.
- **Code:** `STRT-seq_analysis.r:280-356`.
- **Interpretation:** branched tree of EMP→monocyte/macrophage.
- **Reusable?** partial → sc-trajectory (note as **legacy monocle2**; prefer monocle3/slingshot; keep the dispersion-ordering-genes idea).

### M5. Cross-platform / cross-dataset Seurat CCA integration
- **Purpose:** integrate STRT-seq + 10x (and embryonic + external adult datasets) to compare cell types across platforms.
- **How used:** per-batch `NormalizeData` + `FindVariableFeatures(vst, 3000)` → `FindIntegrationAnchors`/`IntegrateData(dims=1:30)`; validate by checking the same DEG markers reproduce across platforms (STRT vs 10x side-by-side dotplots).
- **Code:** `data integration of STRT-seq and 10x.R:60-110`; cross-dataset `STRT-seq_analysis.r:1083-1205`.
- **Interpretation:** platform-robust cell identities; reproducibility across STRT/10x is the validation.
- **Reusable?** partial (dedup) → sc-integration already covers anchors/Harmony/scVI; note the **cross-platform reproducibility check** as a validation pattern only.

### M6. SCENIC binary-regulon consumption (per-cluster regulon activity)
- **Purpose:** summarize TF regulon activity per cluster from a precomputed SCENIC run.
- **How used:** read `binaryRegulonActivity` matrix → drop redundant regulons (Spearman cor > 0.3 with >1 other) → keep regulons active in ≥30% of ≥1 cluster → per-cluster mean activity → pheatmap (ward.D2) with cluster/site annotations.
- **Code:** `STRT-seq_analysis.r:135-168`.
- **Reusable?** yes → sc-grn (SCENIC as the regulon alternative to the house decoupleR route; this is the *consume + filter + per-cluster heatmap* half).

### M7. Plate-based (STRT-seq) normalization & interactive curation
- **Purpose:** correct full-length plate data and manually remove contaminating/doublet cells.
- **How used:** `NormalizeData(scale.factor = 1e5)` (note: 1e5 for full-length STRT/Smart-seq vs 1e4 for 10x); `ScaleData(vars.to.regress = c("nGene","nUMI", batch))`; interactive `DimPlot(..., do.identify=TRUE)` / `CellSelector` to gate out epithelial/mesenchymal/doublet cells on the embedding.
- **Code:** `STRT-seq_analysis.r:11-47`; `data integration ...:` (`CellSelector`).
- **Reusable?** partial → sc-preprocessing (plate-based scale.factor note; interactive gating as a manual-curation option, not the flag-don't-drop default).
- **Anti-pattern to avoid:** the scripts overwrite `@reductions$pca@cell.embeddings` with UMAP coords so `DimPlot(reduction="pca")` shows UMAP. Do **not** copy — store UMAP in its own reduction slot.

## 3. Figure catalog

### F1. Violin + mean-connecting trend across ordered groups
- **Type:** per-gene violin (`scale="width"`) + jittered points + a line segment connecting consecutive group means, across an ordered axis (specification clusters / developmental stages).
- **What it communicates:** the directional trend of a gene along development, with full distribution at each step.
- **Visual style:** `geom_violin(scale="width")` + `geom_point(size=.5)` + `geom_segment` joining `mean[i]→mean[i+1]`; x-axis often blanked when faceting many genes via `ggarrange`.
- **Code:** `STRT-seq_analysis.r:542-591` (head specification), `:619-659` (vs stage).
- **Reusable template:**
  ```r
  # Violin + connecting line of group means across an ORDERED axis. (source: papers/2020-human-macrophage-dev.md F1)
  library(ggplot2)
  ord <- c("Head_Mac0","Head_Mac1","Head_Mac2","Head_Mac3","Head_Mac4")   # ordered groups
  df  <- data.frame(grp = factor(meta$cluster, levels = ord), expr = as.numeric(data[gene, ]))
  m   <- tapply(df$expr, df$grp, mean)
  ggplot(df, aes(grp, expr, fill = grp)) +
    geom_violin(scale = "width") + geom_point(size = .5) +
    geom_segment(aes(x = seq_along(m)[-length(m)], xend = seq_along(m)[-1],
                     y = m[-length(m)], yend = m[-1]), linewidth = 1, inherit.aes = FALSE) +
    theme_classic() + theme(legend.position = "none") + labs(x = NULL, title = gene)
  ```
- **Reusable?** yes → scientific-plotting (ordered-axis violin trend — strong for development/time-course).

### F2. Gene-program pattern: heatmap + faceted trend lines
- **Type:** (a) pheatmap of z-scored stage/cluster means, genes row-ordered by PAM program (navy–white–red); (b) faceted line plot, one facet per program, all genes' trends + a bootstrapped mean ribbon.
- **What it communicates:** discrete temporal/lineage gene programs and their shapes.
- **Visual style:** `pheatmap(color = colorRampPalette(c("navy","white","red"))(100), border_color = NA, annotation_row = program)`; trend `geom_line(alpha low) + facet_wrap(~program, scales="free_y") + stat_summary(fun.data="mean_cl_boot", geom="smooth")`.
- **Code:** `STRT-seq_analysis.r:739-757` (lines), `:984-993` (heatmap).
- **Reusable template:**
  ```r
  # Program trend lines after PAM clustering of stage/cluster-averaged z-scores. (F2)
  ggplot(prog_long) +                                  # cols: gene, program, group(ordered), expression
    geom_line(aes(group, expression, group = gene, color = program), linewidth = .1) +
    facet_wrap(~program, scales = "free_y") +
    stat_summary(aes(group, expression, group = 1), fun.data = "mean_cl_boot",
                 geom = "smooth", color = "black", alpha = .2, linewidth = 1.1) +
    theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ```
- **Reusable?** yes → scientific-plotting (program heatmap + trend-line panel; pairs with M1).

### F3. Stage/subset highlight on UMAP + composition pie
- **Type:** UMAP with one stage/group colored and all other cells greyed; companion pie of that stage's composition with repelled % labels.
- **What it communicates:** where/when a population appears, plus its quantitative makeup at that stage.
- **Visual style:** **grey out others and draw the highlighted cells last** (reorder rows so highlighted plot on top); `scale_color_manual(values = c(pal, others = "#E6E6E6"))`; `ggpubr::ggpie` + `ggrepel::geom_label_repel`.
- **Code:** `STRT-seq_analysis.r:1054-1080` (`plot_single_stage`), `:861-869` (draw-last ordering).
- **Reusable template:**
  ```r
  # Highlight a subset on a UMAP: grey the rest, draw highlighted LAST. (source: papers/2020-human-macrophage-dev.md F3)
  df$grp <- ifelse(df$stage == "CS15", as.character(df$cluster), "others")
  df <- df[order(df$grp == "others", decreasing = TRUE), ]      # 'others' first -> drawn underneath
  ggplot(df, aes(UMAP1, UMAP2, color = grp)) + geom_point() +
    scale_color_manual(values = c(pal_many(n), others = "#E6E6E6")) + theme_bw()
  ```
- **Reusable?** yes → scientific-plotting (highlight-subset idiom; very general).

### F4. Diffusion-map embedding colored by stage / DPT / gene
- **Type:** DC1×DC2 scatter, fill = stage / DPT value / gene; + beeswarm of DC1 (or DPT) by stage.
- **Visual style:** `pch=21` points; **rainbow ramp** `blue→green→yellow→red` for continuous (see §4 — recommend viridis instead); `ggbeeswarm::geom_quasirandom(groupOnX=FALSE)`.
- **Code:** `STRT-seq_analysis.r:884-947`.
- **Reusable?** partial → sc-trajectory (the destiny DC plot; swap rainbow → viridis).

### F5. Circular dendrogram of cell types
- **Type:** hierarchical dendrogram of cell types (ward.D2 on averaged expression of HVGs), leaves colored by cell-type palette, drawn circular.
- **Visual style:** `dendextend::labels_colors()` + `circlize::circlize_dendrogram(dend_track_height=0.7)`.
- **Code:** `STRT-seq_analysis.r:457-464`.
- **Reusable template:**
  ```r
  # Circular cell-type dendrogram. (source: papers/2020-human-macrophage-dev.md F5)
  hc <- hclust(dist(t(avg_expr[hvg, type_order])), method = "ward.D2")
  d  <- dendextend::set(as.dendrogram(hc), "labels_colors", col_type[labels(as.dendrogram(hc))])
  circlize::circlize_dendrogram(d, dend_track_height = 0.7)
  ```
- **Reusable?** yes → scientific-plotting (cell-type similarity dendrogram, linear or circular).

### F6. Averaged-expression marker heatmap (viridis)
- **Type:** pheatmap of z-scored cluster-mean marker expression, clamped to ±1.5, viridis colors.
- **Visual style:** `colorRampPalette(c("#440154","#21908C","#FDE725"))` (viridis); values clamped `[-1.5,1.5]`; `cluster_rows/cols=FALSE` with manual gene/cluster order.
- **Code:** `STRT-seq_analysis.r:118-132`, `:447-454`.
- **Reusable?** partial (dedup) → scientific-plotting heatmaps already cover marker heatmaps; the **clamp-to-±1.5 + viridis** is a style note.

### F7. HSPC-vs-EMP volcano & proportion-over-stage line
- **Type:** two-group volcano (ggrepel top-10) and a cluster-proportion-vs-stage line plot.
- **Code:** `STRT-seq_analysis.r:188-214`.
- **Reusable?** no (dedup) → scientific-plotting volcano.md + composition already cover these.

## 4. Style observations
- **Avoid the rainbow ramp.** The paper leans on a jet-like `blue→green→yellow→red` (and `#E6E6E6→blue→green→yellow→red`) for continuous expression/pseudotime. It's eye-catching but perceptually non-uniform and not colorblind-safe — **prefer viridis** (the paper itself uses `#440154→#21908C→#FDE725` elsewhere). Worth a house rule.
- **Two-tone expression FeaturePlots:** light grey `#E6E6E6` → red/brown (zero expression = light grey) reads cleanly for marker FeaturePlots.
- **Draw-highlighted-last:** reorder rows so the emphasized cells render on top of greyed context — a recurring, generalizable idiom.
- **Connecting-means line over violins/points** to show directional trend across an ordered axis.
- **Clamp z-scores** (±1.5) before heatmaps so a few extreme cells don't wash out the scale.
- **Ordered axes everywhere:** clusters/stages are explicitly releveled to a biological order before every plot.

## 5. Promotion Proposals

| # | target skill | kind | summary | source |
|---|---|---|---|---|
| P1 | scientific-plotting | recipe | **Ordered-axis violin + mean-connecting trend** plot (gene trend across stages/specification steps): violin(scale="width") + points + segment joining consecutive group means. | F1 |
| P2 | scientific-plotting | recipe | **Gene-program pattern panel** — z-scored stage/cluster means → PAM (k via silhouette) → navy–white–red heatmap + faceted trend lines with `mean_cl_boot` ribbon (prep via per-gene ANOVA across ordered groups). | F2, M1 |
| P3 | scientific-plotting (layout-legibility) | recipe | **Highlight-subset-on-UMAP** idiom: grey all others, draw highlighted cells last; + per-group pie with repelled % labels. | F3 |
| P4 | sc-trajectory | method | Add **R-native trajectory options**: destiny `DiffusionMap`+`DPT` diffusion pseudotime, and `princurve` principal-curve lineage fit — alongside the existing PAGA/scanpy route; note monocle2 is legacy (prefer monocle3/slingshot). | M2, M3, F4 |
| P5 | sc-conventions (palettes.md) | requirement | **Continuous-palette rule:** avoid jet/rainbow (blue-green-yellow-red); default to **viridis**; two-tone grey→red OK for marker FeaturePlots; clamp z-scores (±1.5) before heatmaps. | §4, F4, F6 |
| P6 | sc-grn | method | **SCENIC regulon activity** as the alternative to decoupleR: consume binary-regulon output, drop redundant regulons (cor>0.3), keep regulons active in ≥30% of ≥1 cluster, per-cluster mean-activity heatmap. | M6 |
| P7 | sc-preprocessing | requirement | **Plate-based normalization note:** full-length STRT-seq/Smart-seq use `scale.factor = 1e5` + regress nGene/nUMI (vs 1e4 for 10x); record interactive `CellSelector` gating as an optional manual-curation step. Flag the UMAP-into-PCA-slot overwrite as an anti-pattern. | M7 |
