---
paper: "human-pregastrula — Epiblast diversification and blood formation in a human pregastrula"
authors: "ZhenyuXiao lab et al. (s41586-026-10698-y)"
doi: "10.1038/s41586-026-10698-y"
journal: "Nature"
year: 2026
repo: "https://github.com/ZhenyuXiaoLab/Code-of-Human-Early-embryonic-development @ 7c7ff1e0f8edfd9c053337b9862c3ae1ba440ac9"
data_accession: "see paper (spatial transcriptomics + scRNA-seq of human CS6–CS8 embryos)"
distilled: 2026-06-25
tags: [human, embryo, Carnegie-stage, yolk-sac, hematopoiesis, spatial-transcriptomics, scRNA, SCP, scplotter, scop, figure-styling, typography]
---

# human-pregastrula — distillation

## 1. Overview
- **What it is & why notable:** "Epiblast diversification and blood formation in a human pregastrula" (Nature 2026) — a spatial + single-cell study of the human pregastrula (Carnegie stages ~CS6–CS8; epiblast diversification, yolk-sac haematopoiesis, erythroid/megakaryocyte progenitors, spatial slices). Distilled specifically for its **figure styling**: it is a pure figure-code repo (`Figure1–4` + `FigS1–13`, all R) built on the **SCP / scplotter** plotting stack — the same family as our `scop` fork — with consistent, copyable sizing/typography conventions. This is the reference example for "how to get size, font size, theme right."
- **Data:** human; CS6–CS8 embryos; spatial transcriptomics (per-slice, `reduction="spatial"`, up to ~36 slices) + scRNA-seq.
- **Code availability:** **figure code only** (no analysis pipeline, no README/data). 17 scripts, 5,376 lines. Read at snippet fidelity: `Figure2.R` (full), `Figure1.r` (CellDimPlot/integration blocks); styling patterns aggregated by grep across all 17. Stack: SCP + **scplotter** (4 files each), Seurat 4.3, monocle3, slingshot, viridis, ggsci, ggpubr, ggpirate, cowplot, pheatmap.
- **KEY API NOTE:** the repo's `CellDimPlot(obj, group_by=, theme=, palcolor=, raster=)` is the **scplotter/SCP** signature (underscores). Our house fork **scop** uses `group.by=` + `theme_use=` (dots). The *style values* transfer directly; only the argument names differ (see F1 + P2 translation table).

## 2. Methods catalog

### M1. SCP/scplotter plotting stack (vs. our scop fork)
- **Purpose:** publication UMAP/feature/dot plots from Seurat objects.
- **How used:** `CellDimPlot`, `FeatureDimPlot`, `DynamicPlot` (trajectory), `FeatureHeatmap` from SCP/scplotter; plus Seurat `DimPlot`/`DotPlot`/`FeaturePlot` for quick views.
- **Code:** `Figure1.r:161-256`, `:429`.
- **Reusable?** yes (translation) → scientific-plotting: map SCP/scplotter args ↔ scop fork.

### M2. Reference-projection UMAP (query onto reference)
- **Purpose:** place query cells (e.g. CS6 primitive Mk/Ery) onto a reference atlas UMAP to assign identity.
- **How used:** custom `sc_tl_project_umap(ref_data, ref_umap, new_data)` (umap pkg `predict`) → plot reference cells grey + projected cells colored.
- **Code:** `Figure2.R` (project_umap blocks).
- **Reusable?** partial (dedup) → scop already has `ProjectionPlot()`; note as the house equivalent.

### M3. Integration-method benchmark (LISI)
- **Purpose:** compare batch-integration methods (uncorrected / SCT / CCA / Harmony) by iLISI (mixing) and cLISI (cell-type purity).
- **How used:** run each integration, compute iLISI/cLISI, plot; save `5×4`.
- **Code:** `Figure1.r:` (qupici_* blocks, `R2-32_iLISI/cLISI`).
- **Reusable?** partial (dedup) → scop `LISIPlot()`/`BenchmarkPlot()` already cover this; note the 4-method comparison pattern.

### M4. Spatial per-slice mapping & monocle3/slingshot trajectory
- **Purpose:** plot cell types on tissue coordinates per slice; order lineages.
- **How used:** `DimPlot(reduction="spatial", raster=FALSE)` per `slice_num`, patchworked; monocle3 + slingshot for trajectories.
- **Code:** `Figure2.R` (spatial DimPlot), `Figure1.r` / `FigS` (monocle3/slingshot).
- **Reusable?** partial → sc-spatial / sc-trajectory already cover these tools.

## 3. Figure catalog

### F1. Styled UMAP (SCP/scplotter CellDimPlot) — the canonical recipe
- **Type:** UMAP colored by cell type / sample, optionally split by sample.
- **Visual style (the copyable values):** `theme="theme_blank"` (no axes/grid), `pt.size=0.01` (very small for dense embryos), **`raster=FALSE`**, `label=FALSE` with legend on the right, explicit `palcolor` named vector. Single panel saved **8×7**; `split_by="sample"` grid saved **24×21**.
- **Code:** `Figure1.r:161-169, 248-256`.
- **Reusable template (with scop-fork translation):**
  ```r
  # scplotter/SCP signature (as in the paper):
  CellDimPlot(obj, group_by = "celltype", reduction = "umap",
              label = FALSE, palcolor = mycols, pt.size = 0.01,
              theme = "theme_blank", legend.position = "right", raster = FALSE)
  # → our scop fork (rename group_by→group.by, theme→theme_use):
  scop::CellDimPlot(seu, group.by = "celltype", reduction = "umap",
              label = FALSE, palcolor = mycols, pt.size = 0.01,
              theme_use = "theme_blank", legend.position = "right", raster = FALSE)
  # ggsave(..., width = 8, height = 7); split_by="sample" → width = 24, height = 21
  ```
- **Reusable?** yes → scientific-plotting (scop dense-UMAP styling + SCP↔scop arg map).

### F2. Flipped marker dotplot (wide-short)
- **Type:** marker dotplot, genes on y after `coord_flip()`.
- **Visual style:** `DotPlot(..., features=rev(top5)) + coord_flip() + scale_color_gradientn(colours=colorRampPalette(c("#E6E6E6","orange","red","brown"))(100)) + theme(axis.text.x=element_text(angle=60, hjust=1))`. Saved **wide & short: `10 × ~2.0–2.2`** (one short row per gene set).
- **Code:** `Figure2.R` (dotplot), `FigS` (`2ave_*_dot.pdf` saves at `width=10, height=2–2.2`).
- **Reusable?** partial → sizing note (the wide-short aspect for flipped dotplots) → sc-conventions; ramp prefer viridis per house rule.

### F3. Reference-projection UMAP (grey ref + colored query)
- **Type:** UMAP scatter; reference cells `#E6E6E6`, projected query cells colored by predicted type.
- **Visual style:** `geom_point(size=0.8) + theme_bw()`, named palette with `others="#E6E6E6"`.
- **Code:** `Figure2.R`.
- **Reusable?** dedup → scop `ProjectionPlot()`; the grey-ref idiom matches the house highlight convention.

### F4. Cluster-correlation heatmap & GO network
- **Type:** (a) pheatmap of correlation between cluster PCA-mean vectors (`ward.D2`); (b) GO-term igraph network.
- **Visual style:** pheatmap `fontsize_number=8, border_color="grey"`; igraph saved `pdf(width=30, height=30)`, `V(g)$label.cex` scaled (4.5 for GO nodes vs 3.8), bold GO labels.
- **Code:** `Figure2.R` (cor heatmap); `FigS*` (GO graph, `pdf 30×30`).
- **Reusable?** partial → sizing notes (network label.cex + large canvas).

## 4. Style observations (the core of this digest)
Consistent, copyable conventions across all 17 scripts:
- **Global base size:** `theme_classic(base_size = 14)` / `theme_bw()` — 14 pt base everywhere.
- **Axis text:** `element_text(colour = "black", size = 12)`; **titles** `size = 16`, `plot.title = element_text(hjust = 0.5)` (centered).
- **Rotated tick labels:** `angle = 45` or `60` (long names) / `90` (heatmaps), always `hjust = 1`.
- **Point size by density:** `pt.size = 0.01` for dense embryo UMAPs; `0.4–0.5` medium; `1–1.5` sparse. Pair with **`raster = FALSE`** (keep vector) — or rasterize the points layer when n is huge.
- **Cluster label size:** `label.size = 4` on labeled DimPlots; bold (`fontface="bold"`, `family="sans"`).
- **Legend:** default `legend.position = "right"`; `"bottom"` for wide multi-panel; `panel.grid = element_blank()` to de-clutter.
- **Save dimensions by plot type (inches):**
  | plot | width × height |
  |---|---|
  | single UMAP | 8 × 7 |
  | UMAP split by sample (grid) | 24 × 21 |
  | flipped marker dotplot | 10 × ~2.0–2.2 |
  | LISI / metric bar | 5 × 4 |
  | trajectory (slingshot) | 7 × 4.75 |
  | GO igraph network | 30 × 30 |
- **Output:** vector **PDF** as the default device; `viridis`/`ggsci` palettes; `#E6E6E6` light-grey for background/"others".
- **Caveat vs house rules:** the warm `#E6E6E6→orange→red→brown` expression ramp is sequential-OK, but per [[2020-human-macrophage-dev]] P5 prefer **viridis** for continuous and reserve diverging RdBu for signed data.

## 5. Promotion Proposals

| # | target skill | kind | summary | source |
|---|---|---|---|---|
| P1 | sc-conventions (new references/figure-sizing-and-type.md) | requirement | Add a **figure sizing & typography reference**: base_size 14; axis text 12 black; titles 16 centered; rotate 45–60° hjust=1; cluster `label.size≈4` bold; `pt.size` by density (0.01 dense → 1.5 sparse) + `raster` control; legend right/bottom; and a **save-dimensions-by-plot-type table** (UMAP 8×7, split-grid 24×21, flipped dotplot 10×~2.2, metric 5×4, trajectory 7×4.75, network 30×30). | §4, F1, F2, F4 |
| P2 | scientific-plotting (scop-plots.md) | recipe | Add an **SCP/scplotter ↔ scop argument-translation** note so SCP-based published code reuses with our fork: `group_by→group.by`, `theme→theme_use`, plus shared `palcolor`, `split_by`, `pt.size`, `raster`, `label`/`label.size`; include the dense-UMAP styling defaults (theme_blank, pt.size 0.01, raster=FALSE). | F1, M1 |
