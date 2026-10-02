---
name: sc-preprocessing
description: Use when turning CellRanger (or retrieved) counts into a per-sample analysis-ready Seurat object — ambient RNA, doublets, MAD QC; flag, don't drop.
---

# sc-preprocessing (per-sample)

## Overview
CellRanger → ambient removal → Seurat → symbol fix → doublet flagging → QC → (atlas filtering). Runs per library before any merge. **Flags** doublets and low-quality cells (does not delete them), so thresholds can be revisited later. Follows sc-conventions for outputs and logging.

## Entry point
Starting from a public accession instead of your own CellRanger output? Use **sc-dataretrieval** first. True RAW counts enter at step 1; filtered-only matrices use decontX (step 1 alt) or skip ambient removal and flag it.

## Steps
1. **Ambient RNA:**
   - RAW matrix + GPU → **CellBender** on `raw_feature_bc_matrix.h5` → denoised `_filtered.h5`. See references/cellbender.md.
   - Filtered-only / R-only → **decontX** on the filtered matrix. See references/decontx.md.
2. **Load** into Seurat (from CellBender h5 or decontX-corrected counts). For 10x multiome, take `[["Gene Expression"]]`.
3. **Standardize gene symbols** (HGNChelper), collapse duplicates; cross-species ortholog mapping when needed. See references/gene-symbols.md.
4. **Doublets — flag only:** scDblFinder (R) + scrublet (Python), store per-method scores/calls + `doublet_any` and `doublet_consensus`. See references/doublets.md.
5. **QC:** mito/hb/ribo fractions, **MAD-based** thresholds → `qc_pass` flag, pre/post summaries; cell-cycle scoring + regression at normalize/scale. Visualize via scientific-plotting. See references/qc.md.
6. **Atlas gene filtering** (multi-study merges only): keep genes detected in ≥50% of studies + protein-coding; drop tiny samples. See references/atlas-gene-filtering.md.

Log ambient method + params, tool versions, QC cutoffs + rationale, and gene-filtering choices in methods.md.

## Plate-based data (STRT-seq / Smart-seq), not droplet

This pipeline assumes droplet/10x chemistry. For full-length plate-based data, adjust:
(source: https://doi.org/10.1038/s41586-020-2316-7)

- **Skip step 1 (ambient removal) and emptyDrops** — CellBender/decontX/`emptyDrops` assume an empty-droplet
  ambient tail that plate data does not have.
- **Normalization `scale.factor`:** the 1e4 default is a 10x convention; plate libraries have far higher
  counts per cell. Use ~`median(colSums(counts))` (the paper used `1e5`), and regress `nGene`/`nUMI`.
- **Smart-seq (full-length)** benefits from gene-length normalization (TPM-style) for quantitative
  comparisons; **STRT-seq** is UMI-based and does not.
- Doublets: synthetic-doublet classifiers / scDblFinder still apply; QC (MAD) still applies.
- **Anti-pattern (don't copy):** overwriting `seu@reductions$pca@cell.embeddings` with UMAP coords so
  `DimPlot(reduction="pca")` shows UMAP — keep UMAP in its own reduction slot.

## Nuclei (snRNA / single-nucleus): count exon + intron
When you control the count-matrix step for **nuclei** (sci-RNA-seq3, or any re-count from BAM/reads),
count exon and intron features and **sum them** per gene — do not use exon-only.
(source: https://doi.org/10.1126/science.aba7721)

- **Why:** nuclear RNA is mostly unspliced pre-mRNA, so most of a nucleus's signal lives in introns.
  Exon-only counting on nuclei discards the majority of UMIs and cripples sensitivity (Cao et al.
  recovered a UMI/gene count comparable to whole-cell scRNA-seq only because intronic reads were kept).
- **How:** CellRanger `--include-introns` (default in v7+), STARsolo `--soloFeatures Gene GeneFull`
  (`GeneFull` = exon+intron), or velocyto/kallisto spliced+unspliced. If a public **nuclei** matrix was
  built exon-only, treat low UMI/gene yield as a counting artifact, not biology, and re-count if reads exist.
- Whole-cell (droplet) scRNA-seq is dominated by mature mRNA, so exon-only is fine there; this rule is
  specific to nuclei.

## When NOT to use
Integration/merging → sc-integration. Figures → scientific-plotting.
