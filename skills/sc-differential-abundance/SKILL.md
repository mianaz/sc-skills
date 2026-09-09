---
name: sc-differential-abundance
description: Use when testing whether cell-type or neighbourhood proportions shift across conditions (Milo, sccomp, propeller). Not expression DE.
---

# sc-differential-abundance (do proportions shift across conditions?)

## Overview
Tests whether the **relative abundance** of cell types / states changes between conditions — the compositional counterpart to `sc-pseudobulk` (which tests expression). The unit of inference is the **biological replicate (sample), not the cell**. Proportions are **compositional** (they sum to 1, so a rise in one type mechanically lowers others) — a naive t-test/Wilcoxon on raw proportions is statistically invalid; use a method built for count compositions.

## When to use
- "Are there more/fewer <cell type> in disease vs control?" across a cohort of samples
- Subtle or continuous abundance shifts that don't respect cluster boundaries (→ Milo)
- Any composition comparison you intend to put a p-value on

## When NOT to use
- <3 samples per group → no replication; report proportions descriptively, no test
- Expression differences within a cell type → **sc-pseudobulk**
- A stacked-bar/box composition figure with no statistical test → **scientific-plotting**
- Not yet annotated (for cluster-level methods) → **sc-annotation** first (Milo is cluster-free and can run pre-annotation)

## Choose a method
| Method | Level | Use when | Reference |
|---|---|---|---|
| **Milo** (miloR) | neighbourhood (cluster-free) | shifts are subtle/continuous or cluster labels are arbitrary; want to localise *where* on the manifold abundance changes | references/da-methods.md |
| **sccomp** | cluster-level (Bayesian) | want a robust default that models compositional count variability + outliers and gives credible intervals | references/da-methods.md |
| **propeller** (speckle) | cluster-level (frequentist) | few replicates, want a fast classic test (logit/asin transform + limma moderated t/F) | references/da-methods.md |

Decision: cluster-free / localise on manifold → Milo. Cluster-level → sccomp (Bayesian, robust) or propeller (fast, few reps). Reporting more than one is a reasonable robustness check.

## Key caveats
- **Replicates, not cells.** n = samples per group. Thousands of cells from 2 mice is still n=2.
- **Compositionality:** an increase in one type forces apparent decreases elsewhere — interpret shifts jointly, not type-by-type in isolation. Milo/sccomp/propeller account for this; per-type t-tests do not.
- **Annotation/clustering quality** propagates: cluster-level DA is only as good as the labels. Milo sidesteps hard boundaries but still depends on the KNN graph/embedding.
- **Batch confound:** if condition is perfectly confounded with batch/processing day, no method can separate them — check the design first (see sc-pseudobulk on confounding).

## Routes to
- `sc-annotation` — cell-type labels for cluster-level methods
- `sc-pseudobulk` — the expression-shift counterpart (pair the two for "what changed and by how much")
- `scientific-plotting` — composition plots, Milo neighbourhood beeswarm/DA-graph (`references/scop-specialized.md`)
- `sc-conventions` — outputs, source-data export, methods.md logging