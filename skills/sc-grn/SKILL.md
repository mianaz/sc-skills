---
name: sc-grn
description: Use when discovering TF regulons or scoring per-cell TF activity from scRNA-seq (pySCENIC, SCENIC+, AUCell).
---

# Single-Cell Gene Regulatory Networks (pySCENIC / sc-GRN)

## When to use
- Which transcription factors drive each cluster or state
- Per-cell TF activity scores (AUCell regulon matrix)
- De novo TF–target relationships from scRNA-seq
- SCENIC+ (paired RNA + ATAC)

## When NOT to use
- Bulk RNA-seq → WGCNA or `bulk-rnaseq`
- Fewer than ~500 cells → GRNBoost2 is unstable
- GRNBoost2 co-expression only (not the full pipeline) → `arboreto`
- Conserved *gene programs* across species, not TF regulons → `sc-conventions/references/cross-species-tme.md` (cNMF, not pySCENIC)

## Three-stage pipeline
**Stage 1 — Co-expression (GRNBoost2):** TF–target candidates by tree-based regression. Use `arboreto`.

**Stage 2 — Motif validation (cisTarget):** keep pairs whose TF motif is enriched in the target's cis-regulatory regions. Do not skip — co-expression alone is noisy.

**Stage 3 — Activity scoring (AUCell):** per-cell AUC of each regulon. Relative scores; compare within one dataset.

After stage 3, cluster-level on/off summary → [references/regulon-summary.md](references/regulon-summary.md). Motif-DB-free alternative → [references/decoupler-activity.md](references/decoupler-activity.md).

## Key caveats
- cisTarget needs species-specific ranking databases (hg38 or mm10 .feather, ~20 GB); confirm they are on disk before starting.
- SCENIC+ needs paired scATAC; do not run it on RNA-only data.

## Routes to
- `arboreto` — Stage 1 GRNBoost2
- `scanpy` — AnnData that feeds pySCENIC
- `sc-conventions` — cross-species NMF GEP matching
