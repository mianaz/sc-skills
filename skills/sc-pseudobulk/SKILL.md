---
name: sc-pseudobulk
description: Use when comparing scRNA-seq conditions at the sample (not cell) level — pseudobulk DE or sample PCA. Not cluster markers; not cell-type proportions.
---

# sc-pseudobulk (sample-level DE & PCA from single-cell)

## Overview
To compare **conditions** (disease vs control, treated vs untreated, genotype), aggregate each cell's counts up to the **biological replicate** — one pseudobulk profile per sample (optionally per sample × cluster) — then run a bulk DE method (edgeR / DESeq2 / limma-voom). Cells are *measurements*, not replicates; the sample is the unit of statistical inference. This is the field-standard correction to cell-level DE (source: Squair et al. Nat Commun 2021 — cell-level tests treat pseudoreplicates as replicates and inflate false positives; pseudobulk recovers the true DE set).

## When to use
- Differential expression **between conditions**, within a cell type or across all cells
- Sample-level **PCA / MDS / hierarchical clustering** to see how samples group by biology vs batch (your atlas_pseudobulk_pca use case)
- Per-cluster (cell-type-resolved) condition contrasts across a cohort

## When NOT to use
- Marker genes **between clusters in one object** → `Seurat::FindMarkers` / `sc.tl.rank_genes_groups` (that comparison has no biological replicates to aggregate)
- Bulk RNA-seq of whole tissue → `bulk-rnaseq` / `pydeseq2` directly (no aggregation step)
- Fewer than ~3 samples per group → pseudobulk DE has no power; report descriptively, don't fabricate a p-value
- Differential **abundance** (do cell-type proportions shift?) → that is a composition test (Milo/sccomp/propeller), not pseudobulk expression

## Decision tree
1. Comparing conditions with ≥3 samples/group? No → stop (no replication) or use descriptive stats only.
2. Want cell-type-resolved contrasts? Yes → aggregate per **sample × cluster**; No (global) → aggregate per **sample**.
3. Just visualising sample structure, not testing? → aggregate → logCPM → PCA/MDS (skip the DE model).
4. Which DE engine? edgeR-QLF or DESeq2 for counts (default); limma-voom for large cohorts. Method trade-offs → `bulk-rnaseq`; Python DESeq2 → `pydeseq2`.

## Implementation
Aggregation, sample-level PCA, and per-cluster edgeR/DESeq2 contrasts (R-first, Seurat source of truth) → `references/pseudobulk-de.md`.

## Key caveats
- **Sum raw counts, never averages of normalised values** — CPM/logcounts don't re-sum to a valid library.
- Drop sample × cluster cells with too few cells (default: <10 cells) and clusters present in too few samples — thin pseudobulks are unstable.
- Keep sample covariates (batch, sex, donor) in the design matrix; pseudobulk lets you use real bulk covariate modelling.
- `filterByExpr` / independent filtering on the pseudobulk matrix, exactly as in bulk.

## Routes to
- `bulk-rnaseq` — choosing between edgeR / DESeq2 / limma-voom and building the design/contrasts
- `pydeseq2` — running the DESeq2 step in Python
- `scientific-plotting` — PCA/MDS plots and volcano plots from the DE tables (`references/volcano.md`)
- `sc-crispr` — per-perturbation pseudobulk DE is the same aggregation applied to guide identity
- `sc-conventions` — output layout, source-data export, methods.md logging
