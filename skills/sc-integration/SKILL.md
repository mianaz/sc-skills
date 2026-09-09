---
name: sc-integration
description: Use when correcting batch effects or computing integrated embeddings across single-cell samples — Harmony, scVI/scANVI, or scGPT — then carrying results back into the Seurat object.
---

# sc-integration

## Overview
After per-sample preprocessing, merge and integrate. Every path lands embeddings in the Seurat object (the source of truth).

## Choose a path
| Need | Path | Reference |
|---|---|---|
| Fast exploration | Harmony (R) | references/harmony.md |
| Adaptable, downstream-model-compatible; have seed labels → scANVI | scVI/scANVI (Python) | references/scvi-scanvi.md |
| Foundation-model (transformer) embeddings; GPU | scGPT | references/scgpt.md |
| One compartment across many tissues (blood/endo/epi across organs) | per-group downsample + marker-gene-set embedding | references/cross-tissue-embedding.md |
| Across species (human ↔ mouse, …) | shared-ortholog Seurat anchors + k-NN transfer + NNLS no-match gate | references/cross-species.md |

scVI/scANVI and scGPT embeddings are exported from Python and carried back via `CreateDimReducObject` (sc-conventions r-python-bridge). Cluster/UMAP on the integrated reduction.

## QC
Default visual (UMAP by batch vs. biology); scIB/kBET/LISI optional. references/integration-qc.md.

Log method, batch key, latent dims, epochs, rationale in methods.md.

## When NOT to use
Per-sample cleaning → sc-preprocessing. Identity assignment → sc-annotation (scGPT labels are candidates only).
