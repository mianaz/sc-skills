---
name: sc-annotation
description: Use when assigning cell-type identities to single-cell clusters — generating candidate labels via reference-atlas label transfer, pretrained models (SCimilarity, Azimuth), SingleR/celldex, or LLM marker-prompting (GPT-4o/Claude), then confirming every label against canonical (human or mouse) markers with citations using a dotplot, plus AddModuleScore signature scoring (TLS, chemokine, MHC, Ig). Coarse compartments before fine subtypes.
---

# sc-annotation

## Overview
Automation proposes, markers dispose. Candidate labels come from label transfer / pretrained models; **canonical-marker dotplot verification with citations is the mandatory final gate.** Two-level: coarse compartment → fine subtype.

## Flow
1. Cluster on the integrated embedding.
2. Candidate labels (not final) — use one or more, then cross-check agreement:
   - label transfer (references/label-transfer.md)
   - pretrained models (references/pretrained-models.md) — incl. **Garnett** marker-file classifier
     (literature-driven, portable, independent of your clustering)
   - SingleR + celldex reference (references/singler.md)
   - LLM marker-prompting, GPT-4o / Claude (references/llm-annotation.md)
   - **When several references/methods disagree,** quantify agreement with a Cohen's κ
     concordance matrix + spot outlier references → references/label-transfer-concordance.md
3. **MANDATORY:** canonical-marker dotplot verification (references/marker-verification.md) using
   references/markers.md (human) or references/markers-mouse-immune.md (mouse), both cited.
4. Assign identities; log supporting markers + citations + which automated methods agreed in methods.md.
5. **Optional — signature/state scores:** AddModuleScore for published programs (TLS, 12-chemokine,
   MHC-I/II, Ig), or a signature derived from a **bulk** DE table projected onto cells (bulk→sc
   bridge, with pseudobulk reverse-validation). See references/module-scores.md. Compare across
   groups with a test + exact p-value.
6. **Optional — marker-panel QC:** to test whether a small gene panel actually *separates* the
   annotated types (e.g. for a flow panel), measure discriminatory power with a classifier
   benchmarked against ceiling/floor. See references/panel-discriminatory-power.md.
7. **Optional — cluster-specificity QC:** to defend annotation granularity or make keep-vs-merge
   calls, score whether each cluster is separable on the whole transcriptome (SVM cross-validation
   F1 vs a permuted-label null) — clusters near the null floor are over-clustering.
   See references/cluster-specificity.md.

## When NOT to use
Embeddings/clustering inputs → sc-integration. Figures → scientific-plotting.
