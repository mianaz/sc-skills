---
name: sc-crispr
description: Use when analysing pooled CRISPR screens with scRNA-seq readout (Perturb-seq, CROP-seq, Mixscape).
---

# Pooled CRISPR Screens (Perturb-seq / CROP-seq)

## When to use
- User has a Perturb-seq or CROP-seq experiment: scRNA-seq cells each tagged with one or more guide RNAs
- User wants to find which gene knockouts (or activations) drive transcriptional changes
- User wants to cluster perturbations by their transcriptional similarity
- User wants to remove "non-perturbed" cells from guide-assigned cells using Mixscape

## When NOT to use
- Arrayed CRISPR screens (one perturbation per well, bulk RNA-seq readout) → use `bulk-rnaseq` skill instead
- Guide assignment has not been performed → must demultiplex guide identities first (CITE-seq-count or cellranger multi)
- Library coverage < 500 cells per guide → statistical power is insufficient; interpret results cautiously

## Method decision tree
1. **Guide assignment QC:** check UMI counts per guide, fraction of cells with a single vs multiple guides, and fraction of unassigned cells. Flag libraries with >30% unassigned.
2. **Standard scRNA-seq QC and normalisation:** follow `sc-preprocessing` steps on the full object.
3. **Perturbation-aware embedding:** use perturbation-aware PCA (exclude guide-correlated PCs) before UMAP — prevents guide identity from dominating the embedding.
4. **Mixscape (optional but recommended for CRISPRko):** `sc.tl.mixscape` separates cells with effective knockouts from those that escaped editing. Only use for CRISPRko, not CRISPRa/i.
5. **Differential expression per perturbation:** pseudobulk DE comparing each perturbation to non-targeting controls — aggregate per guide/perturbation and run edgeR/DESeq2 via `sc-pseudobulk` (same aggregation, guide identity as the grouping variable).
6. **Perturbation clustering:** use guide-level pseudobulk expression vectors to cluster perturbations with similar transcriptional effects.

## Key caveats
- Guide assignment quality is the most important QC gate — all downstream results are only as good as the assignment
- Non-targeting control cells must be included in the library at sufficient numbers (≥10% of cells recommended)
- Multiple guides per cell (multiplets) confound DE analysis; filter or model them explicitly
- Mixscape requires enough cells per guide to estimate the bimodal KO/non-KO distribution

## Routes to
- `scanpy` — for base scRNA-seq QC, normalisation, and embedding
