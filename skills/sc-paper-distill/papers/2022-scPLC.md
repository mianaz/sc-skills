---
paper: "scPLC — Liver tumour immune microenvironment subtypes and neutrophil heterogeneity"
authors: "Xue R, Zhang Q, Cao Q, et al."
doi: "10.1038/s41586-022-05400-x"
journal: "Nature"
year: 2022
repo: "https://github.com/meta-cancer/scPLC @ 1cd7cac8ebd54088d7bb939c6b51b6715a9d917b"
data_accession: "GSA PRJCA007744 (https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA007744)"
distilled: 2026-06-24
tags: [human, mouse, liver-cancer, HCC, ICC, tumor-immune-microenvironment, neutrophil, TAN, atlas, 10x, Seurat, Scanpy]
---

# scPLC — distillation

## 1. Overview
- **What it is & why notable:** A >1-million-cell liver-cancer atlas (189 samples / 124 patients + 8 mice) that stratifies patients into **five tumour immune microenvironment (TIME) subtypes** (immune activation; suppression by myeloid; suppression by stromal; exclusion; residence) and dissects **tumour-associated neutrophil (TAN)** heterogeneity (CCL4+ recruiting macrophages, PD-L1+ suppressing T cells). Chosen as a first example for its large-atlas engineering (a clean, scriptable preprocessing → dual-toolkit clustering → marker pipeline) and its widely-copied many-cluster UMAP styling.
- **Data:** human + mouse; liver tumours (HCC/ICC and others); 10x scRNA-seq; 189 samples, 124 patients, 8 mice, ~1M cells (Scanpy object observed at 963,276 × 2000 HVG).
- **Code availability:** The repo ships the **analysis pipeline only**, fully:
  - `codes/0-1.read10xCounts.R` — per-sample 10x read-in (DropletUtils).
  - `codes/0-2.prefilter.R` — empty-droplet, QC, doublet removal (scater/scran).
  - `codes/0-3.combine_sclibs.R` — atlas gene-set restriction + merge + normalize.
  - `codes/1a-1`, `1a-2.seurat_cluster*.R` — Seurat clustering (1a-2 adds gene removal).
  - `codes/1b-1..3` — Seurat↔AnnData bridge + `1b-2.scanpy_cluster.py` Scanpy clustering.
  - `codes/2-1.subclus_findmarker.R` — FindAllMarkers + top-N dotplot.
  - `codes/2-2.subclus_marker_plot.Rmd` — lineage-organized marker FeaturePlot panels.
  - **NOT in repo:** every main Nature display figure (the 5 TIME-subtype landscapes, composition heatmaps, chemokine networks, neutrophil trajectories, spatial maps). Those are **caption-only** here; the paper full text is paywalled (no OA), so figure templates for them cannot be faithfully reconstructed without the PDF. See §3 F4.

## 2. Methods catalog

### M1. Empty-droplet removal with a data-driven lower bound
- **Purpose:** drop ambient/empty barcodes per sample before QC.
- **How used:** `DropletUtils::emptyDrops(counts, lower = b)` where `b` = the **second-smallest** per-cell total UMI (`sort(colSums(counts))[2]`) — a data-driven lower bound rather than a fixed number; keep cells at `FDR <= 0.01`.
- **On what data:** each per-sample `SingleCellExperiment` (post `read10xCounts` + `uniquifyFeatureNames`).
- **Tools / packages:** DropletUtils, SingleCellExperiment, BiocParallel.
- **Code:** `codes/0-2.prefilter.R:34-44`.
- **Interpretation:** ambient-aware filtering that adapts the threshold to each library's depth.
- **Reusable?** yes → sc-preprocessing (as an alternative ambient/empty step to the house CellBender/decontX route).

### M2. QC filtering: scater isOutlier + hard caps
- **Purpose:** remove low-quality cells by library size, gene count, and mito fraction.
- **How used:** `isOutlier(sum, log=TRUE, type="lower") | sum > 30000`; `isOutlier(detected, log, lower) | detected < 500 | detected > 6000`; `isOutlier(mito_pct, type="higher") | mito_pct > 50`. Combines adaptive MAD outliers with explicit biological caps. Mito genes from an external `gene_name_type.MT.txt`.
- **Tools / packages:** scater (`perCellQCMetrics`, `addPerCellQC`, `isOutlier`, `plotColData`).
- **Code:** `codes/0-2.prefilter.R:60-75`.
- **Interpretation:** hybrid adaptive+fixed QC; thresholds logged per sample to a `log.file` table.
- **Reusable?** partial → sc-preprocessing (the user's MAD-based QC already covers this; note the hard-cap hybrid + per-sample QC log as a pattern).

### M3. Doublet removal via cluster-level MAD cutoff
- **Purpose:** remove doublets at the cluster level rather than per cell.
- **How used:** `scran::doubletCells()` → per-cell scores; cutoff `= median(log10(scores+1)) + 3*mad(log10(scores+1))`; Louvain-cluster the cells (`buildKNNGraph` k=5 on 50 PCs → `igraph::cluster_louvain`), take the **median score per cluster**, and discard whole clusters whose median exceeds the cutoff.
- **Tools / packages:** scran, scater, igraph, BiocNeighbors.
- **Code:** `codes/0-2.prefilter.R:101-130`.
- **Interpretation:** conservative, cluster-resolution doublet calling (drops doublet-enriched clusters, not scattered cells).
- **Reusable?** yes → sc-preprocessing (alternative to scDblFinder/scrublet; note "flag clusters by median MAD" idea).

### M4. Atlas gene-set restriction before merge
- **Purpose:** standardize the feature space across 189 libraries and drop noise genes.
- **How used:** keep **protein-coding + TR/IG genes only** (`features.prot.tsv` ∪ `gene_name_type.TR_IG.txt`, ~20,234 genes) before merging samples. Separately, `1a-2` removes **MT, HSP, RP, and dissociation-induced** gene sets before HVG selection (external curated lists).
- **Tools / packages:** Seurat, base R.
- **Code:** `codes/0-3.combine_sclibs.R:13-18,40-46`; `codes/1a-2.seurat_cluster_rm.R:30-37,108`.
- **Interpretation:** reduces batch/technical drivers (stress, ribosomal, dissociation) from the variable-gene space — important for a multi-study atlas.
- **Reusable?** yes → sc-preprocessing / sc-integration (corroborates & strengthens the existing protein-coding/multi-study gene-filtering rule; adds the explicit MT/HSP/RP/dissociation removal list).

### M5. Dual Seurat + Scanpy clustering of the same object
- **Purpose:** run the same data through both ecosystems (cross-toolkit reproducibility; pick whichever embedding reads better).
- **How used:** convert with `sceasy::convertFormat` (Seurat→AnnData and back). Both paths use matched params: total-count normalize to **1e4** + log; HVG (Seurat `vst` 2000/1500; Scanpy `n_top_genes`); **regress out percent.mt** on the HVG subset; scale (Scanpy clip at 10); PCA; neighbors (Scanpy `n_neighbors=10`); Louvain (`res` 0.6–1); UMAP. Parameter sweeps are encoded in output filenames (`hvg2000_PC10_res0.6`).
- **Tools / packages:** Seurat, scanpy, sceasy, reticulate, loompy/anndata.
- **Code:** `codes/1a-1.seurat_cluster.R`, `codes/1b-2.scanpy_cluster.py`, `codes/1b-1.scanpy_seurat2anndata.R`.
- **Interpretation:** the same clustering reproduced in two stacks; filename-encoded sweeps make runs self-documenting.
- **Reusable?** partial → sc-conventions (the filename-encodes-params pattern) + r-python-bridge (sceasy round-trip already covered).

### M6. Marker identification + top-N
- **Purpose:** cluster markers for annotation.
- **How used:** `FindAllMarkers(slot='data')` → `group_by(cluster) %>% top_n(5, avg_logFC)`, dedup genes, feed to `DotPlot`.
- **Code:** `codes/2-1.subclus_findmarker.R:55-68`.
- **Reusable?** partial → sc-annotation (standard; the curated marker *sets* in §3 F2 are the higher-value output).

## 3. Figure catalog

### F1. Many-cluster UMAP with a concatenated-qualitative palette
- **Type:** UMAP (DimPlot), categorical color by cluster / sample / type.
- **What it communicates:** global structure across many (>30) clusters and 189 samples without color collisions.
- **Data shown:** UMAP coords; color = cluster id / Sample / Type (Type parsed from sample name via regex `A..._([A-Z_]*)[0-9]?$`).
- **Visual style:** the signature trick — build an **unbounded categorical palette by concatenating every RColorBrewer qualitative palette** (then tripling it), so high-cardinality cluster sets never recycle adjacent colors. `label = TRUE`, `NoLegend()` for the cluster panel, small `pt.size = 0.1` for the dense sample panel.
- **Code:** `codes/2-2.subclus_marker_plot.Rmd` (palette chunk lines ~14-30; plot chunk `DimPlot(..., cols = colorqual)`).
- **Reusable template:**
  ```r
  # High-cardinality categorical palette: concatenate all Brewer qualitative sets.
  # (source: sc-paper-distill/papers/2022-scPLC.md F1)
  library(RColorBrewer)
  qual <- subset(brewer.pal.info, category == "qual")
  colorqual <- unlist(lapply(rownames(qual), function(p)
    brewer.pal(qual[p, "maxcolors"], p)))
  colorqual <- rep(colorqual, 3)            # triple so very large k never runs out
  # sequential stack for ordered/continuous-by-group use:
  seqsets <- c("Blues","Greens","Greys","Oranges","Purples","Reds")
  colorseq <- unlist(lapply(seqsets, function(p)
    brewer.pal(brewer.pal.info[p, "maxcolors"], p)))
  colorseq <- rep(colorseq, 3)
  # Seurat:  DimPlot(seu, group.by = "clusters", cols = colorqual, label = TRUE) + NoLegend()
  ```
- **Interpretation:** lets a single UMAP carry dozens of clusters legibly.
- **Reusable?** yes → sc-conventions/palettes.md (a `pal_many()` high-cardinality categorical palette, alongside the existing Paired ramp).

### F2. Lineage-organized marker FeaturePlot panels + curated liver/TIME marker sets
- **Type:** grid of UMAP FeaturePlots (`ncol = 3`, `label = TRUE`), grouped by immune/stromal/parenchymal lineage.
- **What it communicates:** which clusters express canonical markers, justifying the major cell-type calls.
- **Data shown:** per-gene expression on the UMAP, one small multiple per marker.
- **Visual style:** `fig.width=20, fig.height=12`; literate organization — each lineage block is preceded by a markdown header documenting the marker→cell-type rationale.
- **Code:** `codes/2-2.subclus_marker_plot.Rmd` (FeaturePlot chunks per lineage).
- **Curated liver-TIME canonical markers (directly reusable):**
  - T/NK: `CD3D CD3E CD3G CD8A NKG7 PTPRC`; MAIT/CD4: `CD4 SLC4A10 RORA RORC CXCR4`
  - B/plasma: `MS4A1 CD37 CD79A CD38 MZB1` (+ SSR4)
  - Neutrophil/Mast: `CSF3R S100A8 S100A9 TPSAB1 CPA3 KIT` (+ CXCL8 CXCL2)
  - Macrophage: `CD68 CD14 FCGR3A MARCO C1QC SPP1`
  - Monocyte/DC: `FCN1 CD163 LILRA4 CLEC9A CD1C LAMP3`
  - Endothelial/Fibroblast: `MGP VWF PLVAP CLEC4G COL1A1 ACTA2`
  - Hepatocyte: `GPC3 TTR AFP ALB APOE APOA1`
  - Cholangiocyte: `EPCAM KRT7 PROM1 SOX9 KRT19`
- **Reusable template:**
  ```r
  # Lineage-grouped marker FeaturePlot panel. (source: papers/2022-scPLC.md F2)
  features <- c("CD3D","CD3E","CD3G","CD8A","NKG7","PTPRC")
  Seurat::FeaturePlot(seu, features = features, reduction = "umap",
                      ncol = 3, label = TRUE)
  # house equivalent: scop::FeatureDimPlot(seu, features = features, ncol = 3)
  ```
- **Interpretation:** standard but the **marker sets** are the asset — a vetted liver/TIME panel from a Nature atlas.
- **Reusable?** yes → sc-annotation (liver/TIME canonical marker reference) + scientific-plotting (the panel recipe, noting the scop equivalent).

### F3. Top-N marker DotPlot
- **Type:** DotPlot (mean expression = color, % expressing = size) over top-5 markers per cluster.
- **Visual style:** x-axis labels rotated 90°, `vjust = 0.5`; wide canvas (`width = 20`).
- **Code:** `codes/2-1.subclus_findmarker.R:64-68`.
- **Reusable template:**
  ```r
  # Top-N marker dotplot. (source: papers/2022-scPLC.md F3)
  top5 <- markers |> dplyr::group_by(cluster) |> dplyr::top_n(5, avg_log2FC)
  top5 <- top5[!duplicated(top5$gene), ]
  Seurat::DotPlot(seu, features = top5$gene) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  # NOTE: Seurat ≥4 uses avg_log2FC (this paper's code uses the old avg_logFC).
  ```
- **Reusable?** partial → scientific-plotting (already covered by scop FeatureStatPlot; keep as a Seurat-native fallback + the avg_logFC→avg_log2FC caveat).

### F4. Main Nature display figures — caption-only, NOT in repo
- The five TIME-subtype landscapes, composition/abundance heatmaps, chemokine-network (likely circos/network) panels, neutrophil trajectory, and any spatial maps are **not** in the codebase, and the paper is paywalled (no OA full text retrievable).
- **Status:** `caption-only — no code provided; full text inaccessible`. No reliable template can be reconstructed from the abstract alone.
- **Action if wanted:** user supplies the PDF / figure captions → re-distill into this section (append, don't overwrite).

## 4. Style observations
- **Palette engineering for high cardinality:** rather than a single ramp, concatenate *all* Brewer qualitative sets (and triple) so dozens of clusters get distinct, non-adjacent hues. Mirrors but extends the suite's current `pal_categorical` (Paired ramp). For ordered groups they stack sequential Brewer sets (`colorseq`).
- **Literate marker documentation:** the marker Rmd interleaves markdown headers stating the marker→cell-type logic with the plots — the figure *is* the annotation rationale. Good template for reproducible annotation reports.
- **Self-documenting runs:** every clustering output encodes its params in the filename (`hvg2000_PC10_res0.6`), and preprocessing writes a per-sample QC `log.file` table. A lightweight provenance habit worth adopting.
- **Dense-UMAP hygiene:** small `pt.size` (0.1) and `NoLegend()` on the high-k cluster panel; labels on-data.

## 5. Promotion Proposals

| # | target skill | kind | summary | source |
|---|---|---|---|---|
| P1 | sc-conventions (palettes.md) | palette | Add `pal_many()` — high-cardinality categorical palette by concatenating all Brewer qualitative sets (+ tripling); + sequential stack `colorseq`. Place alongside existing Paired ramp; note when to use which. | F1, §4 |
| P2 | sc-annotation | method | Add a curated **liver / tumour-immune-microenvironment canonical marker reference** (the 8 lineage sets in F2) with citation to Xue et al. Nature 2022. | F2 |
| P3 | scientific-plotting | recipe | Add lineage-grouped **marker FeaturePlot panel** recipe (Seurat `FeaturePlot ncol=3 label`; scop `FeatureDimPlot` equivalent) + Seurat-native top-N DotPlot fallback with the `avg_logFC→avg_log2FC` caveat. | F2, F3 |
| P4 | sc-preprocessing | method | Record an **alternative scran/scater/DropletUtils preprocessing stack** as an option: `emptyDrops` with data-driven lower bound (2nd-smallest UMI, FDR≤0.01) + `doubletCells` cluster-median-MAD cutoff. Mark as alternative to the house CellBender/decontX + scDblFinder route. | M1, M3 |
| P5 | sc-preprocessing / sc-integration | requirement | Strengthen the atlas gene-filtering rule with the explicit removal lists used here: keep protein-coding + TR/IG; drop MT/HSP/RP/dissociation-induced genes **before HVG**. | M4 |
| P6 | sc-conventions | requirement | Adopt **self-documenting run provenance**: encode key params in output filenames (`hvg…_PC…_res…`) and write a per-sample QC log table. (Low priority.) | M5, §4 |
