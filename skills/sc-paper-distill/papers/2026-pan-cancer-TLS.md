---
paper: "Pan-cancer TLS atlas — Pan-cancer spatial atlas of tertiary lymphoid structures"
authors: "Cho KS et al."
doi: "10.1126/science.adz2742"
journal: "Science"
year: 2026
repo: "https://github.com/JialiYue/Pan-Cancer-TLS @ c30154f"
data_accession: "not stated in abstract; paywalled"
distilled: 2026-07-07
tags: [Homo sapiens, pan-cancer, spatial-transcriptomics, Visium, TLS, immune-microenvironment, deconvolution, STRIDE, ssGSEA, survival, ICB-response]
---

# Pan-cancer TLS atlas — distillation

## 1. Overview
- **What it is & why notable:** First pan-cancer spatial atlas of tertiary lymphoid structures (TLSs) across 12 cancer types. Notable for: (1) a computational TLS detection pipeline ("TLS Finder") combining STRIDE deconvolution with B-to-T colocalization scoring on Visium grids, (2) TLS maturation state classification via integration + clustering, (3) an AI framework predicting TLS maturation from H&E, and (4) a maturation-aware composite score that outperforms conventional TLS metrics for patient stratification. The spatial analysis methods are directly reusable for any tissue with immune aggregates.
- **Data:** Homo sapiens, 12 cancer types (BRCA shown in code), 10x Visium spatial transcriptomics, with matched scRNA-seq references for deconvolution. TCGA bulk RNA-seq for survival validation. ICB cohort data from Tiger ICB datasets.
- **Code availability:** The repo ships 7 analysis scripts covering the full pipeline from TLS detection through survival analysis. **Missing:** the AI/deep-learning H&E prediction framework (mentioned in abstract, no code in repo), figure-polishing code, and the ultrahigh-plex single-cell spatial profiling analysis. Plotting is functional but not publication-polished.

## 2. Methods catalog

### M1. B-to-T colocalization scoring on Visium hexagonal grid
- **Purpose:** Quantify spatial co-occurrence of B cells and T cells at each Visium spot to identify candidate TLS locations.
- **How used:** For each spot, compute a colocalization score as `T_prop × B_prop` at the spot itself, plus 0.5× the cross-product with each of its 6 hexagonal neighbors (Visium 55µm layout). Spots passing mean + 0.5×SD for colocalization AND individual T/B proportion thresholds (mean − 0.5×SD) are flagged as TLS candidates. Isolated single spots are removed via hierarchical clustering (`hclust(dist(coor, "maximum"), "single"), h=1`).
- **On what data:** STRIDE deconvolution output (per-spot cell-type fractions) from Visium spatial data.
- **Tools / packages:** Base R, Seurat (for downstream), STRIDE (external deconvolution).
- **Code:** `B_to_T_Colocalization_01.R:40-81` — the hexagonal neighbor lookup + colocalization scoring + thresholding + cluster-based cleanup.
- **Interpretation:** Spots with high B-T colocalization in spatially contiguous clusters represent candidate TLS regions.
- **Reusable?** Yes → sc-spatial. Spatial neighborhood colocalization scoring on Visium hex grids is general-purpose — applicable to any two cell types (e.g., fibroblast–macrophage niches in liver).

### M2. TLS Finder — H&E image segmentation + spot-to-region assignment
- **Purpose:** Map Visium spots to tissue regions using H&E image processing, then assign TLS labels to spots within detected lymphoid aggregates.
- **How used:** Convert H&E to grayscale → multi-Otsu threshold (manual low/high per sample, e.g. 0.2–0.4) → morphological cleanup (remove small objects <500px, fill holes <600px, erosion disk=5) → label connected components → rescale to spot coordinates → expand labels by 50px → assign each spot to its enclosing region.
- **On what data:** H&E images from Visium + tissue_positions_list.csv.
- **Tools / packages:** scikit-image (threshold_multiotsu, morphology, measure, segmentation), matplotlib, pandas.
- **Code:** `TLS_Finder_02.ipynb` cells 6–12. Key parameters are sample-specific (thresholds, object sizes, expansion distance).
- **Interpretation:** Provides a spatial segmentation that anchors molecular TLS calls to histological regions.
- **Reusable?** Yes → sc-spatial. The pattern of H&E segmentation → spot assignment is reusable for any Visium tissue region annotation (tumor vs stroma, necrosis, etc.), though thresholds need per-sample tuning.

### M3. TLS object integration across studies with Harmony
- **Purpose:** Create a unified TLS atlas by integrating TLS-spot pseudo-bulk profiles across multiple BRCA cohorts.
- **How used:** Per TLS: average gene expression across constituent spots → create Seurat object with one "cell" per TLS → merge across cohorts → SCTransform (regress percent.mt) → PCA → Harmony (batch = orig.ident) → UMAP (dims 1:12) → FindNeighbors → FindClusters (resolution 1.4).
- **On what data:** TLS spot-level Visium counts aggregated to TLS-level means.
- **Tools / packages:** Seurat, Harmony, SCTransform.
- **Code:** `Integrate_tls_04.R:357-432` (TLS averaging), `tls_integrated_analysis_05.R:580-608` (integration pipeline).
- **Interpretation:** Harmony-corrected clustering reveals TLS maturation states across cancer types.
- **Reusable?** Partially → sc-integration. The pattern of spot-aggregation → pseudo-bulk → integration is not standard but clever for spatial structure-level analysis. The Seurat+Harmony pipeline itself is standard.

### M4. TLS maturation marker heatmap with ComplexHeatmap
- **Purpose:** Visualize marker gene expression across TLS maturation clusters with biologically meaningful row groupings.
- **How used:** FetchData for curated marker list → compute per-cluster mean expression → Z-score normalize per gene → ComplexHeatmap with row_split by functional category (T cells / B cells / Initiating / Activated), custom color ramp (`colorRamp2(c(-1, 0, 5), c("#118ab2", "#fdffb6", "#e63946"))`), and row annotation sidebar with category colors.
- **On what data:** TLS maturation clusters (4 states).
- **Tools / packages:** ComplexHeatmap, circlize, Seurat.
- **Code:** `tls_integrated_analysis_05.R:622-685`.
- **Interpretation:** Shows progression from naive/initiating to activated TLS states.
- **Reusable?** Yes → scientific-plotting. The ComplexHeatmap template with row_split + category annotation bar is a clean, reusable recipe for any marker-by-cluster heatmap.

### M5. ssGSEA-based ICB response prediction
- **Purpose:** Test whether TLS maturation signatures predict immunotherapy response across ICB cohorts.
- **How used:** Per TLS cluster: extract marker genes (top 100 by avg_logFC) → run GSVA/ssGSEA (method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE) on bulk RNA-seq from ICB cohorts → compare scores between responders (R) vs non-responders (NR) via Wilcoxon test → summarize as median difference (Cor_R) with significance.
- **On what data:** Tiger ICB Dataset bulk RNA-seq across ~25 studies.
- **Tools / packages:** GSVA, Seurat (FindAllMarkers), rstatix (wilcox_test).
- **Code:** `tls_ICB_heatmap_07.R:460-515`.
- **Interpretation:** TLS maturation state 4 (most activated) best predicts ICB response.
- **Reusable?** Partially → could feed a general ssGSEA-based clinical validation recipe, but the ICB-specific framing limits direct reuse for the user's liver atlas work.

### M6. TCGA survival analysis with cluster-correlation approach
- **Purpose:** Validate that TLS maturation signatures correlate with patient survival in TCGA.
- **How used:** Per TLS cluster: compute cluster-average expression in the ST data → correlate each TCGA sample's expression of those genes with the cluster average (Pearson) → use correlation as a continuous score → survival analysis (Cox/KM, implied but not fully shown).
- **On what data:** TCGA RPKM expression + clinical survival data.
- **Tools / packages:** survminer, survival, qpcR.
- **Code:** `tls_survival_analysis_06.R:716-765`.
- **Interpretation:** Higher correlation with activated TLS cluster associates with better survival.
- **Reusable?** Partially → the correlation-as-score approach is an interesting alternative to ssGSEA for mapping spatial signatures to bulk data, but niche.

## 3. Figure catalog

### F1. Spatial B-to-T colocalization spot plot with TLS candidates
- **Type:** Spatial scatter plot (Visium spot layout).
- **What it communicates:** Which spots have high B-T colocalization and which are flagged as TLS candidates.
- **Data shown:** X/Y = Visium spot coordinates; fill = B_T colocalization score (GnBu palette, limits 0–0.3); diamond overlay (pch=5, cyan #00afb9) marks TLS candidate spots.
- **Visual style:** GnBu sequential palette via `colorRampPalette(brewer.pal(9, "GnBu"))`, white point borders (shape=21, stroke=0.5), clean theme (no grid, no axis text), 6×5 inch output.
- **Code:** `B_to_T_Colocalization_01.R:82-100`.
- **Reusable template:**
```r
# Spatial colocalization spot plot with candidate overlay
# (source: sc-paper-distill/papers/2026-pan-cancer-TLS.md F1)
library(ggplot2); library(RColorBrewer)
pal <- colorRampPalette(brewer.pal(9, "GnBu"))
ggplot(spot_df, aes(x = col, y = -row, fill = score)) +
  geom_point(shape = 21, colour = "white", size = 1.5, stroke = 0.5) +
  geom_point(data = candidates, aes(x = col, y = -row),
             pch = 5, fill = NA, size = 1.5, colour = "#00afb9", stroke = 0.3) +
  scale_fill_gradientn(colours = pal(100), limits = c(0, max_val)) +
  labs(x = NULL, y = NULL, fill = "Colocalization") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_blank(), axis.ticks = element_blank())
```
- **Interpretation:** Reader identifies TLS-enriched zones as cyan-marked clusters on a continuous colocalization heatmap.
- **Reusable?** Yes → scientific-plotting. Generic Visium spatial overlay with continuous fill + categorical marker overlay.

### F2. TLS maturation marker heatmap (ComplexHeatmap, row-split)
- **Type:** Heatmap with row splits and annotation bar.
- **What it communicates:** Gene expression patterns distinguishing TLS maturation states.
- **Data shown:** Rows = marker genes (Z-scored mean expression), columns = TLS clusters; row splits by functional category; left annotation bar color-codes category.
- **Visual style:** 3-color gradient (blue #118ab2 → yellow #fdffb6 → red #e63946), row/column gaps (3mm), no column names, row names on right, category colors (#e76f51 T cells, #00CDAC B cells, #f4a261 Initiating, #0077b6 Activated).
- **Code:** `tls_integrated_analysis_05.R:622-685`.
- **Reusable template:**
```r
# ComplexHeatmap with row-split categories and annotation sidebar
# (source: sc-paper-distill/papers/2026-pan-cancer-TLS.md F2)
library(ComplexHeatmap); library(circlize)
col_fun <- colorRamp2(c(-1, 0, 5), c("#118ab2", "#fdffb6", "#e63946"))
category_colors <- c("Cat A" = "#e76f51", "Cat B" = "#00CDAC",
                     "Cat C" = "#f4a261", "Cat D" = "#0077b6")
ha <- HeatmapAnnotation(
  df = data.frame(Category = row_split_factor), which = "row",
  col = list(Category = category_colors))
Heatmap(zscored_mat, col = col_fun,
        cluster_columns = FALSE, cluster_rows = FALSE,
        show_column_names = FALSE, show_row_names = TRUE,
        column_names_side = "top", row_names_side = "right",
        row_split = row_split_factor, column_split = col_split_factor,
        row_gap = unit(3, "mm"), column_gap = unit(3, "mm"),
        left_annotation = ha,
        heatmap_legend_param = list(title = "Z-score",
                                    legend_height = unit(3, "cm")))
```
- **Interpretation:** Clear separation of maturation programs across TLS states.
- **Reusable?** Yes → scientific-plotting. General-purpose ComplexHeatmap recipe with row-split + annotation sidebar.

### F3. ICB response heatmap with significance stars (pheatmap)
- **Type:** Clustered heatmap with significance annotations.
- **What it communicates:** Which TLS maturation clusters predict ICB response across cancer types.
- **Data shown:** Rows = ICB cohorts/cancer types, columns = TLS clusters; fill = median score difference (R vs NR); cell text = significance stars (** p<0.01, * p<0.05).
- **Visual style:** Custom 10-color diverging palette (dark red #67001F → white #ffe8d6 → dark blue #053061), cell size 10×10px, white significance stars, clustered rows/columns, no border.
- **Code:** `tls_ICB_heatmap_07.R:549-565`.
- **Reusable template:**
```r
# pheatmap with significance stars overlay
# (source: sc-paper-distill/papers/2026-pan-cancer-TLS.md F3)
library(pheatmap)
pal <- rev(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7',
             '#ffe8d6','#D1E5F0','#4393C3','#2166AC','#053061'))
# Build significance matrix
pmt <- p_matrix
pmt[p_matrix < 0.01] <- '**'
pmt[p_matrix >= 0.01 & p_matrix < 0.05] <- '*'
pmt[p_matrix >= 0.05] <- ''
pheatmap(score_matrix, scale = "none",
         cluster_rows = TRUE, cluster_cols = TRUE, border = NA,
         display_numbers = pmt, fontsize_number = 12, number_color = "white",
         cellwidth = 10, cellheight = 10, color = pal)
```
- **Interpretation:** ICB response prediction strength varies by TLS maturation state and cancer type.
- **Reusable?** Yes → scientific-plotting. The pheatmap + significance overlay pattern is broadly useful for any score-by-condition matrix with statistical testing.

### F4. TLS feature score spatial plot with cluster labels
- **Type:** Spatial scatter with continuous fill + discrete cluster overlay.
- **What it communicates:** TLS marker gene module score across tissue with TLS region boundaries highlighted.
- **Data shown:** Fill = Seurat AddModuleScore for TLS signature genes (MS4A1, CXCR5, SELL, CD19, LTB, CD79B, CD37, CD79A, TCL1A); overlay = per-TLS cluster labels (colored by Label).
- **Visual style:** Spectral-derived gradient (limits -1 to 1.2), black point borders, TLS spots outlined with cluster-specific colors (pch=21, stroke=1).
- **Code:** `Create_tls_object_03.R:288-307`.
- **Reusable template:** Caption-only adaptation — the code is functional but hardcoded. Template reconstructed in F1 style.
- **Interpretation:** TLS module score validates that detected regions correspond to genuine lymphoid aggregates.
- **Reusable?** Partially → sc-spatial. AddModuleScore + spatial overlay is standard Seurat; the dual-layer (continuous + discrete) is the reusable pattern.

## 4. Style observations
- **Palette family:** Heavy use of RColorBrewer (GnBu sequential, Spectral diverging, RdBu diverging) and ggsci. Custom diverging palette for ICB heatmap.
- **Spatial plots:** Consistent use of shape=21 (filled circle with border), white or black stroke, flipped Y-axis (`y = -row`), clean theme with no gridlines.
- **Heatmaps:** ComplexHeatmap for detailed annotations (row splits, sidebars); pheatmap for simpler score matrices with significance overlays.
- **Color for significance:** White text on colored cells for p-value stars — effective for dark-background heatmaps.
- **Annotation strategy:** Candidate/selected spots shown as diamond overlay (pch=5) rather than subsetting or highlighting, preserving spatial context.
- **TLS marker gene set:** MS4A1, CXCR5, SELL, CD19, LTB, CD79B, CD37, CD79A, TCL1A — used as AddModuleScore feature list. Could be cataloged as a reference gene set.

## 5. Promotion Proposals

| # | target skill | kind | summary | source |
|---|---|---|---|---|
| P1 | sc-spatial | method | Add Visium hexagonal-grid neighbor colocalization scoring recipe: for each spot, compute weighted cross-product of two cell-type proportions with 6 hex neighbors, threshold by mean+k×SD, cluster-cleanup isolated spots | M1 |
| P2 | sc-spatial | method | Add H&E image segmentation → spot assignment pipeline: grayscale → threshold → morphological cleanup → label → rescale → expand → assign spots to regions | M2 |
| P3 | scientific-plotting | recipe | Add ComplexHeatmap row-split recipe with category annotation sidebar (colorRamp2 3-point gradient, row_split factor, HeatmapAnnotation sidebar, row/column gaps) | F2 |
| P4 | scientific-plotting | recipe | Add pheatmap + significance stars overlay recipe (custom diverging palette, display_numbers with ** / * / blank p-value matrix, white number_color) | F3 |
| P5 | scientific-plotting | recipe | Add Visium spatial dual-layer plot recipe (continuous fill + discrete categorical overlay via shape/pch, clean spatial theme, flipped Y) | F1 |
| P6 | sc-conventions | requirement | Establish TLS marker gene reference set: MS4A1, CXCR5, SELL, CD19, LTB, CD79B, CD37, CD79A, TCL1A (for AddModuleScore-based TLS detection in any tissue) | §4 |
