---
name: single-cell
description: Use when starting or orienting any single-cell RNA-seq analysis in an R/Seurat-first workflow. Routes to the stage skill and sc-conventions.
---

# Single-Cell (R/Seurat-first) — Router

## Overview
Opinionated single-cell RNA-seq workflow. **The Seurat object is the source of truth.** Python steps (CellBender, scVI/scANVI, scGPT, scrublet, pretrained annotators) are detours; results are carried back into the Seurat object. Every plot is saved in two versions with its source data; every analytical choice is logged in `methods.md`.

**sc-conventions underlies all stages.** Read it first.

Install the companion `biomed-skill` package for `scientific-plotting` and
`scientific-reproducibility`. These shared skills have one owner and are not
copied into sc-skills. Single-cell figure recipes remain in scientific-plotting.

## Workflow order
0. **sc-dataretrieval** (public accession, not your own CellRanger output) — validate, pick deposited format, load at the right entry point. Discovering datasets before you have an accession → `omic-catalog`.
1. **sc-preprocessing** — CellRanger raw → CellBender → load → doublet flagging → QC. Per-sample. Flag, never silently drop.
2. **sc-integration** — merge, then Harmony, scVI/scANVI, or scGPT (GPU). Embeddings back into Seurat.
3. **sc-annotation** — cluster → candidate labels → mandatory canonical-marker dotplot with citations. Coarse compartment → fine subtype.
4. **scientific-plotting** — house plotting style.

## Routing
- Public accession in (validate, format tier, load) → **sc-dataretrieval**; discover/catalogue first → **omic-catalog**
- Ambient RNA, doublets, QC → **sc-preprocessing**
- Batch correction / latent embeddings / scGPT → **sc-integration**
- Cell-type identities → **sc-annotation**
- Any figure → **scientific-plotting**
- Output paths, palettes, Seurat↔AnnData, methods log → **sc-conventions**
- Trajectory / pseudotime / RNA velocity / lineage → **sc-trajectory**
- Ligand-receptor / CellChat → **sc-cellchat**
- TF regulons / AUCell / pySCENIC → **sc-grn**
- Spatial transcriptomics / Visium / niches → **sc-spatial**
- Perturb-seq / CROP-seq / Mixscape → **sc-crispr**
- scRNA + GWAS / scDRS / target prioritisation → **sc-target**
- Sample-level DE / sample PCA → **sc-pseudobulk**
- Cell-type/neighbourhood proportion shifts → **sc-differential-abundance**
- Learn a recipe from a paper + code (call by name) → **sc-paper-distill**
