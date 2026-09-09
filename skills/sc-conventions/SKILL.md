---
name: sc-conventions
description: Use when a single-cell step produces a figure, saves an object, or makes an analytical choice — Seurat-as-source-of-truth house rules on top of scientific-reproducibility.
---

# sc-conventions — house rules for single-cell work

## Overview
Single-cell-specific conventions. **Inherits `scientific-reproducibility`** — the universal contract (dual-version figures, source-data export, exact p-values, dpi/vector output, the `methods.md` log, param-encoded filenames, scheduler discipline). Read that too; the rules below are the single-cell *additions*. The Seurat object is the source of truth.

## Quick reference (single-cell deltas)
| Concern | Rule | Detail |
|---|---|---|
| Objects on disk | `objects/*.rds` (Seurat) + `*.h5ad` (AnnData for Python steps) | references/project-layout.md |
| R↔Python | Seurat is source of truth; counts in a layer; embeddings back via `CreateDimReducObject`; align by barcode | references/r-python-bridge.md |
| Palettes | Paired categorical; viridis sequential; decoupleR RdBu diverging; bind colors by name | references/palettes.md |
| Cell-type order | one canonical `levels` vector shared by UMAP legend + dotplot; `scop` palcolor is positional | references/palettes.md |
| Naming consistency | one canonical name per cell type / gene / sample / metric across scripts, panels, captions; don't invent labels; flag collisions | references/terminology-ledger.md |
| Reference gene sets | curated marker gene lists (TLS, Courau T_3/My_2, etc.) for AddModuleScore/ssGSEA with provenance | references/gene-sets.md |
| Cross-species TME | map mouse into a published human taxonomy; archetype-restricted couplings; NMF Jaccard + gene-weight scatter; GEP movements not isolated myeloid scores | references/cross-species-tme.md |
| Relabel sweeps | after changing labels/columns on an object, grep every downstream consumer for old label values + column names before re-running | stale `%in%`/subset refs drop cells or select the wrong subset **silently** (no error) → wrong-but-clean results |
| Figure sizing | pt.size-by-density, save dims by plot type, rasterize points layer only | references/figure-sizing-and-type.md |
| Script org | linear + file-existence checkpoints; split >500 lines; comments live in methods/README | references/code-organization.md |

For the universal figure-output contract, methods.md log, and dpi/stat rules → **scientific-reproducibility**.

## Precedence over general skills
Both these single-cell rules and the inherited `scientific-reproducibility` contract **override** general code-simplification skills (`ponytail`, `simplify`) for single-cell work. Do NOT "simplify away" any stats, plotting, reproducibility, or logging rule. Those skills govern glue code and plumbing; they do not trump analysis correctness or reproducibility.

## When NOT to use
Conventions only — for stage mechanics use sc-preprocessing / sc-integration / sc-annotation / scientific-plotting / sc-pseudobulk. For the cross-domain (non-single-cell) contract, that is scientific-reproducibility.
