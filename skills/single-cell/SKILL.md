---
name: single-cell
description: Use when starting or orienting a single-cell RNA-seq analysis in an R/Seurat-first workflow — covers the core pipeline (CellBender preprocessing, doublet flagging, Harmony or scVI/scANVI integration, marker-verified annotation) and routes to the stage-specific skill and the shared conventions.
---

# Single-Cell (R/Seurat-first) — Router

## Overview
Opinionated single-cell RNA-seq workflow. **The Seurat object is the source of truth.** Python steps (CellBender, scVI/scANVI, scrublet, pretrained annotators) are run as detours and their results are always carried back into the Seurat object. Every plot is saved in two versions with its source data; every analytical choice is logged in `methods.md`.

This is the **core pipeline backbone** — preprocessing → integration → annotation, plus shared conventions, the reproducibility contract, and figures — with **sc-paper-distill** as the growth engine that extends the suite from the literature. Downstream analyses (trajectory, gene regulatory networks, cell–cell communication, spatial, CRISPR screens, differential abundance, pseudobulk DE) are out of scope for this backbone; add them via the paper-distill → promotion loop.

## Workflow order
1. **sc-preprocessing** — CellRanger raw → CellBender → load → doublet flagging (scDblFinder + scrublet) → QC. Per-sample. Flag, never silently drop.
2. **sc-integration** — merge, then Harmony (fast, R) or scVI/scANVI (adaptable, Python → carried back to Seurat).
3. **sc-annotation** — cluster → candidate labels (label transfer / pretrained) → MANDATORY canonical-marker dotplot verification with citations. Two-level: coarse compartment → fine subtype.
4. **scientific-plotting** — R-first plotting style (scop, tidyplots, ggplot2+cowplot, ComplexHeatmap).

**sc-conventions underlies all stages** — project layout, palettes, dual-version + source-data + stats + sizing/dpi rules, R↔Python bridge, methods-log discipline. Read it first. It layers on top of **scientific-reproducibility**, the cross-domain output contract.

## Dependencies
Each skill carries its core method inline. Outside names that appear are **software packages** (Seurat, scanpy, edgeR, DESeq2, scVI, …) — install the ones a step actually uses; they are libraries, not skills.

## Routing
- Loading raw counts, removing ambient RNA, calling doublets, QC → **sc-preprocessing**
- Batch correction / latent embeddings → **sc-integration**
- Assigning cell-type identities → **sc-annotation**
- Making any figure → **scientific-plotting**
- Output paths, palettes, saving rules, methods log, Seurat↔AnnData conversion → **sc-conventions**
- The universal output / reproducibility contract for any data analysis → **scientific-reproducibility**
- Learning a method/figure from a paper+code and folding it back into the suite → **sc-paper-distill**
