---
paper: "FCA-liver — Decoding human fetal liver haematopoiesis"
authors: "Popescu, Botting, Stephenson et al."
doi: "10.1038/s41586-019-1652-y"
journal: "Nature"
year: 2019
repo: "https://github.com/haniffalab/FCA_liver @ 9abb77b8eddc2dc7c42d722fbc3d45cb417cb66d"
data_accession: "ArrayExpress (developmentcellatlas.ncl.ac.uk portal); confirm E-MTAB id before citing"
distilled: 2026-06-24
tags: [human, fetal-liver, haematopoiesis, development, atlas, 10x, Smart-seq2, FDG, PAGA/AGA, diffusion-map, trajectory]
---

# FCA-liver — distillation

## 1. Overview
- **What it is & why notable:** The Haniffa-lab human fetal liver haematopoiesis atlas — a landmark developmental-haematopoiesis map. Chosen for **method novelty + iconic visuals**: a custom force-directed-graph (FDG/ForceAtlas2) embedding and animation of the differentiation topology, AGA (abstract graph abstraction — the *precursor to scanpy PAGA*), and a clean, reusable R plotting framework (`dr.plot` + separate indexed legends, diagonalized marker heatmap/spotplot). Unlike most repos, this one **ships the figure-generating code**.
- **Data:** human; fetal liver (plus skin/kidney/yolk-sac comparisons in the paper); 10x + Smart-seq2; ~140k cells; scRNA-seq + flow/CITE.
- **Code availability:** unusually complete — a 22-step numbered `pipelines/` workflow + a `tools/` library. Read at snippet fidelity: `tools/bunddle_utils.R` (plotting + DR + classifier utils), `pipelines/09_gene_heatmap_and_spotplot/gene_heatmap_and_spotplot.R`, and `README.md` (documents every step). Described from README (not line-read): FDG animation (`pipelines/14_*`), pseudotime (`pipelines/13_*`), random-forest discriminatory power (`pipelines/12_*`), SVM classifiers (`pipelines/15/16`), batch correction (`pipelines/19` Harmony, `06` ComBat), clustering comparison (`pipelines/18`), AGA/cell comparison (`pipelines/22`). Seurat **v2.3.4** era — APIs (`@data`, `@dr`, `SetAllIdent`) are legacy; templates below are modernized.

## 2. Methods catalog

### M1. Force-directed graph (FDG / ForceAtlas2) embedding for differentiation topology
- **Purpose:** lay out continuous differentiation trajectories better than t-SNE/UMAP for branching haematopoiesis; also animate it.
- **How used:** ForceAtlas2 (custom iterative `fa2`) run on a PCA + shared-nearest-neighbor (SNN) graph; ~600 iterations; ~2 h for 100k cells. `runFDG(pca.df, snn, iterations, ...)` writes PCA+SNN to disk, calls `make_fdg.py`, returns X/Y coords. Animation via OpenCV over FA2 iterations.
- **Tools / packages:** custom `fa2` (modified forceatlas2), `Matrix::writeMM`; OpenCV for video.
- **Code:** `tools/bunddle_utils.R:5-30` (runFDG); `pipelines/14_fdg_animation_write_input/*` (animation, README-level).
- **Interpretation:** the iconic fetal-liver FDG shows HSC/MPP → erythroid/myeloid/lymphoid branches as a connected topology.
- **Reusable?** yes → sc-trajectory (FDG as a topology-preserving layout; modern equivalent below).

### M2. AGA (abstract graph abstraction) — precursor to PAGA
- **Purpose:** abstract the single-cell graph into a cell-type/cluster graph with connectivity weights (which populations are continuous vs. distinct).
- **How used:** scanpy's early AGA implementation; "AGA score" plots between cell types; hard to read with >10 types (per README).
- **Code:** `pipelines/11_plot_dr/plot_dr.R`, `pipelines/22_cell_comparison/cell_comparison.R` (README-level).
- **Interpretation:** quantifies trajectory connectivity between annotated populations.
- **Reusable?** yes (as guidance) → sc-trajectory: **AGA is the precursor to `scanpy.tl.paga`; use PAGA today** (often to initialize FDG/`draw_graph`).

### M3. Gene-set discriminatory power via random forest
- **Purpose:** quantify how well a chosen gene set separates cell types — a *QC metric for marker/flow panels*, not a production classifier. Chosen for conceptual similarity to flow-cytometry gating.
- **How used:** random forest on a user gene set, downsampled to 5000 cells/type; outputs classification report + confusion matrix.
- **Tools / packages:** scikit-learn RandomForest (Python).
- **Code:** `pipelines/12_gene_discriminatory_power_analysis/*` (README-level).
- **Interpretation:** high accuracy ⇒ the gene panel is sufficient to discriminate those populations (e.g. for designing a flow panel).
- **Reusable?** yes → sc-annotation / sc-target: validate a marker panel's discriminatory power.

### M4. Diffusion-map + monocle pseudotime
- **Purpose:** order cells along a single developmental lineage and find pseudotime-varying genes.
- **How used:** `destiny` diffusion map (a DC used as pseudotime) refined with `monocle` 2.6.4. Strong precondition: input must be a *single true lineage* (doublets/outliers removed first), or diffusion maps give nonsense.
- **Tools / packages:** destiny 2.6.2, monocle 2.6.4.
- **Code:** `pipelines/13_pseudotime/pseudotime.R`, `heatmap_plot.R` (README-level).
- **Reusable?** partial → sc-trajectory (the single-lineage precondition is the reusable lesson; the suite's modern stack is scVelo/monocle3/slingshot).

### M5. SVM classifiers — cell-type & doublet (synthetic-doublet training)
- **Purpose:** label transfer (cell type) and doublet flagging via learned models.
- **How used:** SVM on top-20 DEGs/type (PCA features); doublet SVM trained on *synthetic* doublets vs. real cells. `Apply_Classifier_On_Seurat_Object()` pads missing feature genes with zeros so a classifier can run on a new object. README stresses classifiers are **tissue-specific** (don't cross tissues).
- **Code:** `tools/bunddle_utils.R:78-150` (apply); `pipelines/15_train_classifier`, `16_train_doublets_classifier` (README-level).
- **Reusable?** low (dedup) → sc-preprocessing/sc-annotation already cover doublets (scrublet/scDblFinder) and label transfer; note the synthetic-doublet idea + tissue-specificity caveat only.

## 3. Figure catalog

### F1. Indexed DR plot — numbered population labels at medians (`dr.plot`)
- **Type:** 2D DR scatter (UMAP/t-SNE/FDG/diffusion), colored by cell label.
- **What it communicates:** many populations on one embedding without label collisions — each population gets a **numbered gray disc at its median**, decoded by a separate legend.
- **Data shown:** DR1/DR2; color = population; numeric index = population (1..n).
- **Visual style:** palette = `sample(colorRampPalette(brewer.pal(12,"Paired"))(n))` with `set.seed` (shuffled Paired ramp so adjacent clusters differ); `pt.size = 0.2`; label discs `pch=21, colour="gray", alpha=.5` at `aggregate(median)` positions; optional 2nd categorical var as point **shape** (circle vs triangle).
- **Code:** `tools/bunddle_utils.R` (`dr.plot`, ~L210-285; `dr.plot.numerical` for continuous).
- **Reusable template:**
  ```r
  # Indexed DR plot: number each population at its median; decode via separate legend.
  # (source: sc-paper-distill/papers/2019-FCA-liver.md F1 — Popescu et al. Nature 2019)
  library(ggplot2); library(dplyr)
  df <- data.frame(DR1 = emb[,1], DR2 = emb[,2], lab = factor(labels))
  med <- df |> group_by(lab) |> summarise(DR1 = median(DR1), DR2 = median(DR2))
  med$idx <- seq_len(nrow(med))
  ggplot(df, aes(DR1, DR2, color = lab)) +
    geom_point(size = 0.2) +
    geom_point(data = med, color = "gray50", fill = "gray", shape = 21, size = 6, alpha = .5) +
    geom_text(data = med, aes(label = idx), color = "black", size = 3) +
    scale_color_manual(values = pal_many(nlevels(df$lab))) +  # sc-conventions pal_many
    theme_classic() + theme(legend.position = "none")          # legend rendered separately, F2
  ```
- **Interpretation:** compact, panel-friendly cell-type map.
- **Reusable?** yes → scientific-plotting (indexed-label DR plot recipe; scop `CellDimPlot` is the house default, this is the numbered-index variant for very many populations).

### F2. Separate indexed legend PDF (`plot.indexed.legend`)
- **Type:** standalone legend (grid of numbered color swatches + text), `theme_void`.
- **What it communicates:** the index→population key, rendered as its **own vector PDF** so figure panels can be assembled/re-laid-out freely.
- **Visual style:** `ncols` grid; swatch `geom_point(size=symbol.size)`; index number + label text annotations; nothing else.
- **Code:** `tools/bunddle_utils.R` (`plot.indexed.legend`, ~L190-210).
- **Reusable template:**
  ```r
  # Standalone legend as its own vector PDF (decouples legend from plot for panel assembly).
  # (source: papers/2019-FCA-liver.md F2)
  plot_indexed_legend <- function(labels, colors, ncols = 2, symbol = 8, text = 4) {
    n <- length(labels); nr <- ceiling(n / ncols)
    d <- data.frame(X = rep(seq_len(ncols), each = nr)[seq_len(n)],
                    Y = rep(rev(seq_len(nr)), times = ncols)[seq_len(n)],
                    idx = seq_len(n), lab = labels, col = colors)
    ggplot(d, aes(X, Y)) + geom_point(size = symbol, colour = d$col) +
      geom_text(aes(label = idx), size = text) +
      geom_text(aes(x = X + .1, label = lab), hjust = 0, size = text) +
      theme_void() + xlim(0.5, ncols + 1)
  }
  ```
- **Reusable?** yes → sc-conventions (publication figure-assembly convention: export legends as separate vector PDFs).

### F3. Mean-expression heatmap + spotplot with **diagonalization**
- **Type:** marker heatmap (ggplot `geom_tile`) + companion dot/"spot" plot (size+color = expression).
- **What it communicates:** which markers mark which cell types, with genes **reordered so high expression runs down the diagonal** — far more readable than input order.
- **Data shown:** rows = cell types, cols = genes; fill/size = mean expression per type.
- **Visual style:** heatmap `scale_fill_gradient(low="lightblue", high="darkred")`, black tile borders, x labels 45°; spotplot `scale_color_gradient(low="lightsteelblue1", high="darkred")` + matched `scale_size_continuous`; legends merged via `guide_legend()`.
- **Code:** `pipelines/09_gene_heatmap_and_spotplot/gene_heatmap_and_spotplot.R` (heatmap ~L120-140; spotplot ~L150-175; diagonalize ~L100-115).
- **Reusable template (the diagonalization algorithm — the real asset):**
  ```r
  # Diagonalize a celltype x gene mean-expression matrix: order genes by weighted
  # center-of-mass so high values form a diagonal. (source: papers/2019-FCA-liver.md F3)
  # mat: rows = cell types (in desired order), cols = genes; values = mean expression.
  weighted_center <- function(v) sum(seq_along(v) * v) / sum(v)
  gene_order <- order(apply(mat, 2, weighted_center))
  mat <- mat[, gene_order]
  # alt ordering: hclust(dist(t(mat)), method = "ward.D")$order
  # then melt -> geom_tile (heatmap) or geom_point(aes(size, color)) (spotplot).
  ```
- **Interpretation:** a self-organizing marker panel; the diagonal makes type-specific genes obvious.
- **Reusable?** yes → scientific-plotting (add diagonalization to heatmaps.md / complexheatmap.md — applies to any marker heatmap/dotplot).

### F4. Force-directed graph embedding + 2D animation
- **Type:** FDG scatter (X/Y from ForceAtlas2) + MP4 animation over FA2 iterations.
- **What it communicates:** the branching topology of differentiation as a physical layout; the animation shows cells settling into lineage arms.
- **Visual style:** same `dr.plot` indexed-label styling on FDG coords; precomputed `label_colours.csv` color key; legend PNG.
- **Code:** `tools/bunddle_utils.R:5-30` (runFDG) + `pipelines/14_fdg_animation_write_input/*` (README-level). Example movie: developmentcellatlas.ncl.ac.uk/datasets/liver_fdg_movie/.
- **Reusable template (modern equivalent):**
  ```python
  # Modern FDG: PAGA-initialized force-directed layout in scanpy. (papers/2019-FCA-liver.md F4)
  import scanpy as sc
  sc.pp.neighbors(adata, n_neighbors=15)
  sc.tl.paga(adata, groups="cell_type")          # AGA's modern successor (M2)
  sc.pl.paga(adata, plot=False)
  sc.tl.draw_graph(adata, init_pos="paga")        # ForceAtlas2 layout, PAGA-seeded
  sc.pl.draw_graph(adata, color="cell_type")
  ```
- **Reusable?** yes → sc-trajectory (FDG/`draw_graph` topology layout, PAGA-initialized).

### F5. Pseudotime expression curves
- **Type:** gene expression vs pseudotime scatter + loess smoothing curve (normalized and raw).
- **What it communicates:** how marker genes rise/fall along a developmental lineage.
- **Code:** `pipelines/13_pseudotime/pseudotime.R`, `heatmap_plot.R` (README-level; template reconstructed).
- **Reusable template (reconstructed):**
  ```r
  # Gene-vs-pseudotime curve. (source: papers/2019-FCA-liver.md F5, reconstructed)
  df <- data.frame(pt = pseudotime, expr = as.numeric(expr_gene))
  ggplot(df, aes(pt, expr)) +
    geom_point(size = .2, alpha = .3) +
    geom_smooth(method = "loess", se = TRUE, color = "darkred") +
    labs(x = "pseudotime", y = "expression") + theme_classic()
  ```
- **Reusable?** partial → sc-trajectory/scientific-plotting (the suite's scvelo/trajectory plots cover most of this).

### F6. Interactive WebGL viewers + AGA graph (mostly superseded)
- **Type:** standalone HTML 2D/3D WebGL embeddings with hover labels; interactive AGA graph.
- **Status:** functional but **deprecated by the authors** in favor of web portals (Fast Portals). Heavy embedded-data HTML, internal-sharing only.
- **Code:** `tools/interactive_2D_viewer/*`, `interactive_3D_viewer/*`, `gene_expression_viewer_apps/*`.
- **Reusable?** no (low value; modern portals — cellxgene, ShinyCell, Vitessce — supersede). Recorded for completeness.

## 4. Style observations
- **Index-and-decode labeling:** instead of crowding cell-type names onto the embedding, number populations at their medians and ship the key as a *separate* legend — keeps dense embeddings clean and panels re-arrangeable. Pairs with the suite's `pal_many()` for high cardinality.
- **Legends as first-class separate artifacts:** every DR plot saves a standalone legend PDF — a deliberate figure-assembly workflow (mix/match panels in Illustrator without baked-in legends).
- **Diagonalization as default for marker matrices:** order genes by weighted center-of-mass so the heatmap/dotplot self-organizes into a diagonal — adopt for any marker heatmap.
- **Shuffled + seeded categorical palette:** `sample(colorRampPalette(Paired)(n))` with a fixed seed — shuffling separates neighbors; the seed keeps it reproducible. (Superseded by `pal_many()` from [[2022-scPLC]], but same intent.)
- **Naming hygiene:** alphanumeric-only cell-type/cluster names — special chars (`/ \ @`) break the downstream tooling. A small but real reproducibility rule.
- **Two-tone expression ramp:** light→`darkred` (lightblue / lightsteelblue1 → darkred) for mean-expression heat/spot plots.

## 5. Promotion Proposals

| # | target skill | kind | summary | source |
|---|---|---|---|---|
| P1 | scientific-plotting (heatmaps.md / complexheatmap.md) | recipe | Add **gene diagonalization** for marker heatmaps/dotplots: order genes by weighted center-of-mass (`Σ i·vᵢ / Σ vᵢ`) so expression runs down the diagonal (+ hclust ward.D alternative). | F3 |
| P2 | scientific-plotting | recipe | Add the **indexed DR-plot** recipe: number populations at their medians (gray disc + index) for very-high-cardinality embeddings, as the numbered-index variant of scop `CellDimPlot`. | F1 |
| P3 | sc-conventions (plot-rules.md) | requirement | Adopt **separate vector legend PDFs** as a figure-assembly convention (decouple legend from plot for panel layout); + **alphanumeric-only label naming** hygiene. | F2, §4 |
| P4 | sc-trajectory | method | Add **PAGA-initialized force-directed graph** (`scanpy draw_graph(init_pos="paga")`) as a topology-preserving layout for differentiation; note **AGA = PAGA precursor** and the diffusion-map single-lineage precondition. | M1, M2, M4, F4 |
| P5 | sc-annotation / sc-target | method | Add **gene-set discriminatory-power QC**: random forest (downsampled per type) to test whether a marker/flow panel separates cell types (flow-gating analogy). | M3 |
