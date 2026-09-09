---
paper: "HuMu TME — Differential assembly of mouse and human tumor microenvironments"
authors: "Courau T et al."
doi: "10.1038/s41590-026-02505-7"
journal: "Nature Immunology"
year: 2026
repo: "not available — no public analysis/figure repo; HTO helper scripts only at https://github.com/UCSF-DSCOLAB/aarao_scripts/; browser at https://quipi.org/app/quipi_humu"
data_accession: "GSE310560 (mouse scRNA-seq); human ImmunoProfiler as Combes et al. Cell 2022; 2,543 mouse bulk RNA-seq from SRA/recount3"
distilled: 2026-09-04
tags: [Mus musculus, Homo sapiens, pan-cancer, TME, CyTOF, scRNA-seq, cNMF, chemokine, immune-archetypes, cross-species]
---

# HuMu TME — distillation

## 1. Overview
- **What it is & why notable:** Courau, Combes, Fragiadakis, and Krummel systematically ask whether the mouse models that dominate immuno-oncology actually assemble the same tumor microenvironment (TME) as human tumors. Across 15 commonly used models (CyTOF in all; scRNA-seq in nine), most murine TMEs map onto a minority of human tumors: immune-desert, macrophage-rich archetypes (~17% of ImmunoProfiler patients). Chemokine networks, T–myeloid frequency coupling, and most gene-expression-program (GEP) “movements” diverge — but a core set of T-cell and myeloid GEPs is conserved, and one coordinated axis (T-cell cytotoxicity × myeloid IFN response) is conserved and prognostic. Notable for *how* they compare species (composition embedding, chemokine compartment scaling, cNMF + Jaccard + gene-weight scatter, GEP–GEP movements), not only for the biological punchline. Directly relevant to any cross-species cancer atlas that currently frames mouse models as generally representative.
- **Data:** Mouse: 15 models (B16F10, MC38, CT26, LLC, 4T1, RENCA, KPC, MMTV-PyMT, and others), BALB/c and C57BL/6, mostly day-14 transplanted tumors (~300–500 mm³); CyTOF *n* = 109 tumors; scRNA-seq in 9 models after lymphoid/myeloid FACS (59,551 cells in the merged UMAP). Human: UCSF ImmunoProfiler (flow + bulk RNA-seq of sorted T, Treg, myeloid, tumor, stroma; composition *n* ≈ 224; imaging *n* = 85) plus a 13-indication TCGA subset (*n* = 4,341) that does **not** include prostate. External: 2,543 mouse bulk RNA-seq samples / 178 SRA studies classified into Thorsson immune subtypes.
- **Code availability:** **not available** for analysis or figures. Public pieces: GEO `GSE310560` (processed Seurat RDS noted on seqout), interactive browser https://quipi.org/app/quipi_humu, HTO demultiplexing helper https://github.com/UCSF-DSCOLAB/aarao_scripts/, source-data XLSX for plotted numbers, Supplementary Tables 1–3 (sample metadata, GEP top-50 genes, antibodies). Figure code, cNMF scripts, and cosine-similarity pipeline were not released. Evidence level: `full-text` (OA) for methods/results/captions; `not available` for repo.

## 2. Methods catalog

### M1. Ten-feature composition embedding and cosine similarity to human TME archetypes
- **Purpose:** Put mouse models and human tumors in one composition space so “does this model look like human?” is a quantitative mapping, not a qualitative impression.
- **How used:** Build a 10-feature immune matrix (PanCan features from Combes et al. 2022: CD45+, T, Treg, CD4 Tconv, CD8, non-granulocytic myeloid, monocytes, macrophages, cDC1, cDC2, stroma — as used in Fig. 1c). Z-score each feature across the pooled 224 human + 109 mouse samples. Aggregate humans by ImmunoProfiler archetype and mice by tumor line. Embed with UMAP; hierarchically cluster; compute cosine similarity (row-normalize to unit length, then `S = X X^T`; cosine distance = 1 − S).
- **On what data:** Mouse CyTOF frequencies vs human flow-inferred frequencies (human flow parameters inferred from gene scores as in ref. 7). Imaging of T:myeloid ratios used as orthogonal confirmation (QuPath + StarDist, random-forest classifier).
- **Tools / packages:** UMAP; cosine similarity in base linear algebra; ggplot2/ComplexHeatmap implied; QuPath 78 + StarDist 79 for imaging.
- **Code:** not provided. Methods paragraph “Comparison of human and mouse TME diversity” is the spec.
- **Interpretation:** Nearly all models (except RENCA) cluster with immune-desert, macrophage-rich human archetypes that are only ~17% of ImmunoProfiler patients (0% of HCC to 46.5% of gynecologic tumors). Mouse variance is myeloid-driven; human variance is T-cell-driven. T:myeloid imbalance survives orthotopic vs subcutaneous KPC, “dirty” housing, aging, and high-fat diet.
- **Reusable?** yes → sc-differential-abundance / sc-conventions. Any cross-species atlas should map mouse samples into a published human TME taxonomy rather than only showing within-atlas proportion bars.

### M2. ImmuneSubtypeClassifier of 2,543 mouse bulk tumors onto Thorsson C1–C6
- **Purpose:** Test whether the desert bias is an artifact of 15 lab models or a property of published mouse tumor RNA-seq at large.
- **How used:** MetaSRA query (RNA-seq, mouse, tissue, cancer DOIDs) → manual filter excluding PDX/sorted/organoid/non-tumor → 2,543 samples / 178 studies. recount3 counts, 75th-percentile + log2 normalize, Babelgene mouse→human orthologs, ImmuneSubtypeClassifier (BestCall). Also: z-score Thorsson signature genes in TCGA, take PC1 gene weights, apply those weights to mouse, rank medians 1–6 per species (Extended Data Fig. 1h).
- **On what data:** SRA mouse tumor bulk RNA-seq vs Thorsson et al. Immunity 2018 human subtypes (TCGA).
- **Tools / packages:** ontologyIndex, MetaSRA, recount3, Babelgene, ImmuneSubtypeClassifier (Gibbs bioRxiv 2020).
- **Code:** not provided.
- **Interpretation:** >60% of mouse samples call as lymphocyte-depleted (C4); ~7% as immune-rich (C2/C6). Low T-cell infiltration is the default published mouse tumor, not a quirk of B16/MC38.
- **Reusable?** yes → sc-target / bulk deconvolution recipes. Useful as an external check on whether a GEMM cohort is typical-mouse or an outlier (RENCA-like).

### M3. Compartment-scaled chemokine ligand/receptor comparison
- **Purpose:** Ask whether the cell-recruitment grammar that builds TMEs is the same across species.
- **How used:** Human: TMM-CPM, log, then *z*-score each chemokine/receptor **across compartments** (T, Treg, myeloid, tumor, stroma) so the heatmap shows which compartment owns the gene, not absolute TPM. Mouse: scRNA-seq pseudobulk mean per sample per compartment, then the same across-compartment scaling. Orthologs via BioMart (April 2018 archive). Protein follow-up for CCR5 by flow. Human scRNA-seq of HNSC used to show CXCL13 is usually T-cell-biased but occasionally stromal.
- **On what data:** ImmunoProfiler bulk-sorted compartments vs mouse scRNA-seq (including FACS-contaminating tumor/stroma kept because ratios were stable). Mouse “stroma” defined by scoring the 20-gene human stroma signature from Combes 2022 on non-immune cells; Fibroblasts 1/2 and Myofibroblasts 1 called stroma.
- **Tools / packages:** edgeR-style TMM (human); Seurat pseudobulk (mouse); BioMart; flow cytometry.
- **Code:** not provided.
- **Interpretation:** Conserved: CCR4, CCR7, CXCR4, CXCR6 and their ligands. Divergent: CCR2/CCR5 transcripts (and CCR5 protein) lower on mouse T cells than human T cells. CXCR3/CXCR5 receptors conserved on T cells, but ligands relocate — murine CXCL9/CXCL10 and CXCL13 are stromal (fibroblast), whereas human CXCL9/10 are myeloid-biased and CXCL13 is T/Treg-biased. Direct implication: CXCL13+ T-cell / TLS / ICB-response circuits studied in patients rarely exist in standard mouse models.
- **Reusable?** yes → scientific-plotting + a chemokine method recipe. The *across-compartment z-score, species side-by-side* design is the reusable unit, not a generic heatmap of chemokine genes.

### M4. Sample-level cell-frequency correlations with mouse-like-archetype restriction
- **Purpose:** Conserved *abundances* can hide non-conserved *couplings*. Test which cell–cell density relationships survive the species jump.
- **How used:** Pearson correlation matrices of immune feature frequencies (BH-adjusted). Plot human matrix hierarchically clustered; plot mouse matrix in the **human column order**. For discordant pairs, recompute the human correlation after restricting to ImmunoProfiler archetypes that resemble mice (ID CD4/CD8 mac). Conserved examples: immune infiltrate vs tumor Ki-67 (inverse); CD8 frequency vs exhaustion (PD-1+CTLA-4+); cDC1/cDC2 vs CD4/CD8. Discordant globally: CD4 Tconv–Treg; T–myeloid; CD8–macrophage — restored when humans are censored to mouse-like deserts.
- **On what data:** Same 10-feature matrices as M1, plus CyTOF Ki-67 and exhaustion gates (Extended Data Fig. 4).
- **Tools / packages:** Pearson + Benjamini–Hochberg; ggplot2 scatter with `lm` + s.e. ribbon.
- **Code:** not provided.
- **Interpretation:** Mouse T–macrophage coupling is a property of desert TMEs, not of human cancer in general. Drug studies that perturb macrophages or CD8 cells in MC38/B16 should be read as modeling that subset. Guardrail: do not treat a mouse correlation as a human universal without the archetype-restricted sensitivity analysis.
- **Reusable?** yes → sc-differential-abundance. Cheap, high-leverage panel for any paired human/mouse atlas.

### M5. Consensus NMF of T-cell and myeloid GEPs, Jaccard matching, gene-weight scatter, external MP validation
- **Purpose:** Move from “which cell types are present” to “which transcriptional programs are the same molecule.”
- **How used (human bulk):** TPM from sorted conventional T cells (Treg excluded) and non-granulocytic myeloid; EHK ≥ 8; log2; top 5,000 genes by MAD. Rank sweep 3–30, 50 NMF runs, cophenetic correlation coefficient (CCC) on consensus connectivity; local CCC maxima as candidate ranks; 10 extra runs concatenated, *k*-means on H10, local outlier factor (40% cutoff), median consensus W/H (DECIPHER-seq-like). 9 stable T GEPs, 14 myeloid GEPs.
- **How used (mouse scRNA):** DECIPHER-seq `iNMF_ksweep()` on conventional T and non-granulocytic myeloid, `batch = FALSE`, *k* = 2–40 × 20 replicates, sample as grouping variable. 25 T GEPs, 23 myeloid GEPs (two apparent myeloid duplicates noted).
- **Matching:** Jaccard on top 20 genes (T) or top 50 (myeloid); Fisher exact test. Gene-weight scatter of human vs mouse loadings for matched pairs (dashed line = 40th highest contributor). External Jaccard vs Gavish et al. 2023 meta-programs and Miller et al. 2025 glioma myeloid/T GEPs; also vs Hu et al. 2023 wound-healing myeloid GEPs. Mouse GEP usage heatmaps: scale H across cells, average by subtype.
- **On what data:** ImmunoProfiler bulk T/myeloid vs HuMu mouse scRNA-seq T/myeloid (granulocytes excluded). Authors flag bulk-vs-scRNA as likely **underestimating** conservation (Extended Data Fig. 7b).
- **Tools / packages:** sklearn NMF; DECIPHER-seq 57; cophenetic CCC (Brunet 2004); local outlier factor; Jaccard + Fisher.
- **Code:** not provided.
- **Interpretation:** Conserved T: T_3 cytotoxicity (PRF1, LAG3, GZMB, NKG7) and T_9 CD4 regulation (TSC22D3, JUNB, RGS1, IL7R, CD69). Conserved myeloid: My_1 inflammatory (IL1A, IL1B, NLRP3), My_2 IFN response (IFIT2, IFIT3, ISG15, CXCL10), My_8 lipid-metabolism TAMs, My_11 LYVE1 TAMs (LYVE1, FOLR2, CD163, MRC1, SELENOP). These match independent human meta-programs. GEP identity is a better cross-study/cross-species currency than cluster labels.
- **Reusable?** yes → sc-grn / existing NMF pipelines. Steal: (i) split T vs myeloid before NMF, (ii) Jaccard *and* gene-weight scatter, (iii) Fisher test, (iv) external MP Jaccard, (v) drop contamination GEPs (their T_2/T_6/My_10). Do not copy the bulk-vs-scRNA design if both species are scRNA.

### M6. Intercellular GEP “movements” and 2×2 survival
- **Purpose:** Treat coordination between a T-cell program and a myeloid program as the translatable unit, not either program alone.
- **How used:** Pearson correlation of T-GEP vs myeloid-GEP enrichment across samples (BH). The only movement clearly conserved in mice is T_3 cytotoxicity ↔ My_2 IFN response (visible in IR CD8 mac / IS CD8 human archetypes and in MC38, CT26, B16-F10). Bin patients high/low by median (ImmunoProfiler) or top/bottom 50% (TCGA text) / top/bottom 30% (TCGA methods — **internal inconsistency**; use one and report it). Four-group KM + log-rank. TCGA GEP score = percentile-normalize top 20 contributor genes, then average.
- **On what data:** ImmunoProfiler GEP usage; TCGA 13-type subset (BLAD, CRC, GBM, GYN, HNSC, KID, HEP, LUAD, PDAC, SARC, MEL — **no PRAD**).
- **Tools / packages:** Pearson + BH; survival / log-rank; ggplot2 KM.
- **Code:** not provided.
- **Interpretation:** T_3Hi My_2Hi is the most favorable OS combination in ImmunoProfiler (trend) and TCGA (clearer). T_9 CD4-regulation ↔ My_1 inflammatory is **not** conserved in mice; T_9Hi My_1Lo trends better in TCGA; mouse-like (T_9Lo, macrophage-rich desert) archetypes are among the worst survivors. Framing: mice are most relevant to the human patients with currently poor outcomes.
- **Reusable?** yes → survival + NMF. Highest-leverage figure upgrade for a myeloid-focused cancer atlas: stop scoring myeloid programs in isolation; test the coordinated T×myeloid axis.

## 3. Figure catalog

### F1. Fig 1c–d — species violins of TME composition and relative variance
- **Type:** grouped violins (1c) + variance scatter (1d).
- **What it communicates:** Mice are not a noisy version of humans; they are systematically T-poor, macrophage-rich, and their diversity lives in myeloid composition.
- **Data shown:** 1c: frequency of named compartments (out of Live, T, or myeloid as labeled) in human (blue, *n* = 224) vs mouse (orange, *n* = 109); Bonferroni-corrected two-sided *t*-tests with exact *P*. 1d: relative variance of each parameter, human vs mouse, diagonal = equal variance, points colored by compartment.
- **Visual style:** two-color species palette (human blue / mouse orange) used consistently for the whole paper; Tukey box overlay implied in later boxplots; asterisks with exact *P* in legend; no rainbow cell-type coloring at this stage.
- **Code:** caption-only — no code provided; template reconstructed in house style.
- **Reusable template:**
```r
# Species-split composition violins
# (source: sc-paper-distill/papers/2026-humu-tme.md F1)
library(ggplot2)
ggplot(freq_df, aes(species, frequency, fill = species)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white") +
  facet_wrap(~ feature, scales = "free_y") +
  scale_fill_manual(values = c(human = "#4C78A8", mouse = "#F58518")) +
  labs(x = NULL, y = "Frequency") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", strip.background = element_blank())
```
- **Interpretation:** The translational default mouse TME is the human desert, not the human mean.
- **Reusable?** yes → scientific-plotting. Species-split violin/box with a locked two-color palette.

### F2. Fig 2a–b — mouse models embedded in human archetype space
- **Type:** UMAP of composition features + hierarchically clustered cosine-similarity heatmap.
- **What it communicates:** Mapping, not averaging: which model ≈ which human archetype.
- **Data shown:** 10 z-scored feature frequencies; human samples as circles colored by archetype (*n* = 224), mouse lines as triangles (*n* = 12 after averaging); heatmap of cosine similarity mouse line × human archetype.
- **Visual style:** human = circles, mouse = triangles (shape encodes species so color can encode archetype); heatmap with species color bars (human blue, mouse orange); hierarchical clustering on both axes for 2a-right, cosine matrix for 2b.
- **Code:** caption-only — no code provided; template reconstructed.
- **Reusable template:**
```r
# Cosine similarity of mouse models to human TME archetypes
# (source: sc-paper-distill/papers/2026-humu-tme.md F2)
library(ComplexHeatmap); library(circlize)
X <- scale(feature_mat)           # samples x 10 features, pooled z-score
Xn <- X / sqrt(rowSums(X^2))      # cosine = unit-row product
S <- tcrossprod(Xn)
mouse_ids <- rownames(S)[species == "mouse"]
human_ids <- rownames(S)[species == "human"]
Heatmap(
  S[mouse_ids, human_ids],
  name = "cosine",
  col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
  rect_gp = gpar(col = "white", lwd = 0.4),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2"
)
```
- **Interpretation:** RENCA is the exception (immune-rich, CD4/Treg/DC-high, T_3/My_2-low, T_9-high). Everyone else sits on desert/macrophage archetypes.
- **Reusable?** yes → scientific-plotting. Shape-encodes-species UMAP + cosine heatmap is the figure we should copy for GEMM-to-patient mapping.

### F3. Fig 3a,d,e — chemokine receptor/ligand compartment heatmaps
- **Type:** side-by-side scaled heatmaps + exemplar boxplots.
- **What it communicates:** Receptors can be conserved while ligand *source* is not — assembly logic diverges even when the receptor toolkit looks similar.
- **Data shown:** rows = chemokine receptors (3a) or selected ligands (3d); columns = human vs mouse T cells (3a) or compartments (3d); values = expression scaled across compartments within species. 3e: CCL22, CXCL9, CXCL13 TPM/pseudobulk by compartment.
- **Visual style:** human|mouse column blocks; hierarchical clustering of the full matrix in Extended Data Fig. 3; boxplots with Tukey whiskers, one point per sample.
- **Code:** caption-only — no code provided; template reconstructed.
- **Reusable template:**
```r
# Across-compartment z-score, species side-by-side chemokine heatmap
# (source: sc-paper-distill/papers/2026-humu-tme.md F3)
library(ComplexHeatmap); library(circlize)
scale_within_species <- function(mat) {
  t(scale(t(mat)))  # genes x compartments, z across compartments
}
H <- scale_within_species(human_cpmt)  # e.g. CXCL* x (T, Treg, My, Tumor, Str)
M <- scale_within_species(mouse_cpmt)
Heatmap(
  cbind(H, M),
  name = "z",
  col = colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426")),
  column_split = factor(c(rep("Human", ncol(H)), rep("Mouse", ncol(M)))),
  cluster_columns = FALSE,
  rect_gp = gpar(col = "white", lwd = 0.3)
)
```
- **Interpretation:** Do not use mouse CXCL13+ stroma as a model of human CXCL13+ T cells / TLS.
- **Reusable?** yes → scientific-plotting. Column-split species chemokine heatmap.

### F4. Fig 4a–d — paired correlation matrices and conserved vs discordant scatters
- **Type:** two Pearson heatmaps (human clustered; mouse in human order) + side-by-side scatterplots.
- **What it communicates:** Conservation is a property of *relationships*, and some relationships only exist in the human subset that looks like mouse.
- **Data shown:** cell-type frequencies; inset scatters after restricting humans to ID CD4/CD8 mac archetypes; `lm` line + s.e. ribbon; points colored by archetype (human) or tumor line (mouse).
- **Visual style:** identical axis limits and point size across species columns; conserved pairs on the left, discordant on the right (caption grouping); BH stars on the matrix.
- **Code:** caption-only — no code provided; template reconstructed.
- **Reusable template:**
```r
# Human-ordered correlation matrix pair
# (source: sc-paper-distill/papers/2026-humu-tme.md F4)
library(corrplot)
Rh <- cor(human_freq, use = "pairwise.complete.obs")
ord <- hclust(as.dist(1 - Rh))$order
Rm <- cor(mouse_freq, use = "pairwise.complete.obs")
corrplot(Rh[ord, ord], method = "color", tl.col = "black", tl.cex = 0.7)
corrplot(Rm[ord, ord], method = "color", tl.col = "black", tl.cex = 0.7)
```
- **Interpretation:** CD8–macrophage coupling in GEMMs should be reported as desert-TME biology, then tested in human PCa stratified by immune-rich vs desert.
- **Reusable?** yes → scientific-plotting.

### F5. Fig 5b–g — Jaccard GEP similarity + gene-weight scatters
- **Type:** Jaccard heatmaps (Fisher stars) + human-vs-mouse gene-loading scatterplots.
- **What it communicates:** A matched GEP is not just overlapping top genes; the *weights* of shared drivers should also agree. Species-private drivers sit off-diagonal.
- **Data shown:** Jaccard of top 20 (T) / top 50 (myeloid) genes; scatter of gene weights for T_3, T_9, My_2, My_1; bold = shared high contributors; gray = species-private; dashed lines cut the top 40 contributors.
- **Visual style:** highlighted matched pairs in bold colored text on the Jaccard heatmap; two-panel gene-weight scatters with equal aspect; no density contour (points only).
- **Code:** caption-only — no code provided; template reconstructed.
- **Reusable template:**
```r
# Gene-weight scatter for a matched human/mouse NMF pair
# (source: sc-paper-distill/papers/2026-humu-tme.md F5)
library(ggplot2)
cutoff <- 40
ggplot(w, aes(human_weight, mouse_weight)) +
  geom_point(aes(color = rank_min <= cutoff), size = 1.2, alpha = 0.85) +
  geom_hline(yintercept = sort(w$mouse_weight, decreasing = TRUE)[cutoff],
             linetype = 2, color = "grey50") +
  geom_vline(xintercept = sort(w$human_weight, decreasing = TRUE)[cutoff],
             linetype = 2, color = "grey50") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey70")) +
  coord_equal() +
  labs(x = "Human gene weight", y = "Mouse gene weight", color = "Top 40") +
  theme_classic(base_size = 11)
```
- **Interpretation:** Jaccard-only matching (especially at a 0.05 threshold) over-calls conservation; the scatter is the quality control.
- **Reusable?** yes → scientific-plotting + NMF method. This is the panel our current cross-species NMF figure is missing.

### F6. Fig 6a–f — GEP-movement heatmap, T_3×My_2 scatter, 2×2 Kaplan–Meier
- **Type:** T-vs-myeloid GEP Pearson heatmaps; scatter of two GEP scores; 2×2 KM.
- **What it communicates:** The translatable object is a *movement* (two programs rising together), and it parses survival better than either program alone.
- **Data shown:** 6a human vs mouse GEP–GEP correlation; 6b T_3 vs My_2 scatter (human colored by archetype, mouse by line); 6c median split; 6d ImmunoProfiler OS; 6e–f same logic in TCGA with top-20-gene scores.
- **Visual style:** two-column species heatmaps; scatter with `lm` + s.e.; four KM curves (T_3HiMy_2Hi vs the other three combinations), not a single high/low split.
- **Code:** caption-only — no code provided; template reconstructed.
- **Reusable template:**
```r
# 2x2 GEP combination KM
# (source: sc-paper-distill/papers/2026-humu-tme.md F6)
library(survival); library(survminer)
clin$T3 <- ifelse(clin$T3_score >= median(clin$T3_score), "Hi", "Lo")
clin$My2 <- ifelse(clin$My2_score >= median(clin$My2_score), "Hi", "Lo")
clin$axis <- factor(paste0("T3", clin$T3, "_My2", clin$My2),
                    levels = c("T3Hi_My2Hi", "T3Hi_My2Lo", "T3Lo_My2Hi", "T3Lo_My2Lo"))
fit <- survfit(Surv(time, event) ~ axis, data = clin)
ggsurvplot(fit, pval = TRUE, risk.table = TRUE,
           palette = c("#1B9E77", "#7570B3", "#D95F02", "#666666"))
```
- **Interpretation:** Favorable immunity is interferon-licensed myeloid cells co-occurring with cytotoxic T cells, not myeloid infiltration per se.
- **Reusable?** yes → scientific-plotting. Direct upgrade to an MDSC-only survival figure.

## 4. Style observations
- Locked **two-color species system** (human blue, mouse orange) across violins, UMAPs, heatmaps, and scatters. Archetype or tumor-line color is secondary and only used when species is already encoded by shape or facet.
- **Paired species columns** are the layout idiom: human left, mouse right, shared axis limits, shared feature order (mouse heatmaps forced into human dendrogram order).
- Heatmaps use white cell borders, hierarchical clustering where the claim is unsupervised, and **suppressed clustering** where the claim is “same order, different pattern.”
- Statistics are visible: exact *P* in Fig. 1 captions, Bonferroni or BH named, Tukey box definition written out, two-sided tests declared. Survival uses log-rank, not just a Cox table.
- Scatters always carry `lm` + s.e. ribbon; they never rely on a correlation coefficient in the caption alone.
- Resource framing: a public browser (quipi) is treated as part of the paper, not an afterthought.

## 5. Promotion Proposals

| # | target skill | kind | summary | source |
|---|---|---|---|---|
| P1 | scientific-plotting | recipe | Paired-species composition violins + locked human-blue / mouse-orange palette | F1 |
| P2 | scientific-plotting | recipe | Cosine-similarity heatmap of mouse models (or GEMM genotypes) vs human TME archetypes, with shape-encodes-species UMAP | F2 |
| P3 | scientific-plotting | recipe | Across-compartment z-scored chemokine ligand/receptor heatmap, column-split by species | F3 |
| P4 | scientific-plotting | recipe | Human-ordered dual correlation matrices + conserved/discordant scatter pair | F4 |
| P5 | scientific-plotting | recipe | NMF gene-weight scatter (human vs mouse loadings, top-40 dashed cutoffs) | F5 |
| P6 | scientific-plotting | recipe | 2×2 Kaplan–Meier for a coordinated T×myeloid GEP axis | F6 |
| P7 | sc-conventions | requirement | Cross-species TME papers must map mouse samples into a published human taxonomy (ImmunoProfiler archetypes and/or Thorsson C1–C6) rather than only reporting within-atlas proportions | M1, M2, §4 |
| P8 | sc-conventions | requirement | Conserved cell-type frequency *correlations* require an archetype-restricted sensitivity analysis (recompute in the human subset that resembles mouse) | M4 |
| P9 | sc-grn (or NMF method notes) | method | Split T vs myeloid before cNMF; match with Jaccard **and** gene-weight scatter; Fisher-test overlaps; Jaccard-validate against Gavish/Courau meta-programs; drop contamination factors | M5 |
| P10 | sc-target / survival notes | method | Score coordinated GEP “movements” (T cytotoxicity × myeloid IFN) and stratify survival in 2×2, not as single signatures | M6 |

`kind` ∈ { recipe, requirement, method, palette, fix }.

## 6. Directed distill tables (when requested)

### 6.1 Paper info table

| title | authors | journal | year | doi | main findings summary |
|---|---|---|---:|---|---|
| Differential assembly of mouse and human tumor microenvironments | Courau T et al. | Nature Immunology | 2026 | 10.1038/s41590-026-02505-7 | Most of 15 common mouse tumor models recapitulate immune-desert, macrophage-rich human TMEs (~17% of patients); chemokine ligand sources and many T–myeloid couplings diverge; conserved T-cell cytotoxicity and myeloid IFN GEPs form a prognostic “movement.” Prostate is absent from their TCGA comparator. |

Persisted in `references/paper_catalog.csv`.

### 6.2 Marker/signature curation table

| marker_or_signature | marker_type | experiment | species | context | claim_use | evidence_source | citation |
|---|---|---|---|---|---|---|---|
| T_3 T-cell cytotoxicity (PRF1, LAG3, GZMB, NKG7) | signature | bulk RNA-seq (human sorted T) / scRNA-seq (mouse T) | mouse\|human | conserved T GEP | identity of cytotoxic program; one pole of the conserved movement | Fig. 5c, Extended Data Fig. 5c, Supplementary Table 2 | 10.1038/s41590-026-02505-7 |
| T_9 CD4 regulation (TSC22D3, JUNB, RGS1, IL7R, CD69) | signature | bulk RNA-seq / scRNA-seq | mouse\|human | conserved T GEP | CD4 regulatory circuit; not globally coordinated with My_1 in mice | Fig. 5d, Extended Data Fig. 7e–h | 10.1038/s41590-026-02505-7 |
| My_1 inflammatory (IL1A, IL1B, NLRP3) | signature | bulk RNA-seq / scRNA-seq | mouse\|human | conserved myeloid GEP | IL-1 inflammatory myeloid state; poor-outcome axis when high with low T_9 | Fig. 5g, Extended Data Fig. 5d | 10.1038/s41590-026-02505-7 |
| My_2 IFN response (IFIT2, IFIT3, ISG15, CXCL10) | signature | bulk RNA-seq / scRNA-seq | mouse\|human | conserved myeloid GEP | IFN-stimulated myeloid state; pairs with T_3 for favorable OS | Fig. 5f, Fig. 6 | 10.1038/s41590-026-02505-7 |
| My_8 lipid metabolism TAM | signature | bulk RNA-seq / scRNA-seq | mouse\|human | TAM-associated GEP | correlates with TAM density; list not fully enumerated in main text | Fig. 5e, Extended Data Fig. 6b,j | 10.1038/s41590-026-02505-7 |
| My_11 LYVE1 TAMs (LYVE1, FOLR2, CD163, MRC1, SELENOP) | signature | bulk RNA-seq / scRNA-seq | mouse\|human | resident-like TAM GEP | conserved LYVE1 TAM program; also similar to wound-healing GEP | Extended Data Fig. 6k, Extended Data Fig. 7d | 10.1038/s41590-026-02505-7 |
| CXCL13 | gene/protein | bulk + scRNA-seq | mouse\|human | chemokine ligand source | T/Treg-biased in human TME, stroma-biased in mouse TME | Fig. 3d–f | 10.1038/s41590-026-02505-7 |
| CXCL9 / CXCL10 | gene/protein | bulk + scRNA-seq | mouse\|human | CXCR3 ligands | myeloid-biased in human, stromal in mouse | Fig. 3d–e | 10.1038/s41590-026-02505-7 |
| CCR2 / CCR5 | gene/protein | RNA + flow (CCR5 protein) | mouse\|human | T-cell chemokine receptors | reduced on mouse vs human T cells | Fig. 3a–c | 10.1038/s41590-026-02505-7 |
| CXCR4 / CXCR6 / CCR4 / CCR7 | gene/protein | bulk + scRNA-seq | mouse\|human | T-cell chemokine receptors | conserved expression pattern across species | Fig. 3a, Discussion | 10.1038/s41590-026-02505-7 |
| T:myeloid imaging panel (CD3, CD4, CD8, CD163, HLA-DR, XCR1, EpCAM human; CD3e, CD11b, IA/IE mouse) | modeling marker set | 7-plex IF / IF | mouse\|human | T vs myeloid ratio | orthogonal confirmation that human tumors are T-dominated, mouse tumors myeloid-dominated | Fig. 1e–f, Methods | 10.1038/s41590-026-02505-7 |

Then appended to `references/omic_marker_db.csv`. Full top-50 lists live in Supplementary Table 2 (not ingested here; `list not fully enumerated in manuscript text` for My_8).

## 7. Evidence level tags

| item | evidence |
|---|---|
| M1–M6 methods parameters | `full-text` (Methods + ED Figs) |
| F1–F6 visual style | `captions-only` plus in-text figure descriptions (raster figures not re-derived from PDF) |
| Analysis / figure code | `not available` |
| GEP top-50 full lists | Supplementary Table 2 named, not parsed (`methods-only` / table-not-ingested) |
| Interactive browser / GEO | `full-text` Data availability |
| Prostate absence from TCGA 13-type list | `full-text` Methods “TCGA cohort and survival analyses” |
