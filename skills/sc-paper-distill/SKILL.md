---
name: sc-paper-distill
description: Distill a single-cell, spatial or omics paper and its code into a sourced methods digest and reusable analysis proposals.
disable-model-invocation: true
---
# Single-cell Paper Distillation

Read the paper, relevant supplements, figures and available code. Reuse supplied
local assets first; retrieve missing material when accessible. Record title,
authors, year, DOI and the code commit. State which sources were actually read.

Use [references/digest-schema.md](references/digest-schema.md) for the digest and
[references/extraction-tips.md](references/extraction-tips.md) for source retrieval
and code-to-figure mapping. Distinguish published implementation from reconstructed
examples and biological findings from proposed extensions.

Save the digest in the user's research output directory, not the installed skill:
`papers/<year>-<shortname>.md`. Record methods and figures with stable IDs, exact
source locators, parameters, independent units, controls and interpretation. Add
a dated revision when updating a previous digest.

When requested, create `paper_catalog.csv` and `omic_marker_db.csv` beside the
output digests. These are user-generated outputs, not bundled databases. Follow
[references/omic_marker_db_schema.md](references/omic_marker_db_schema.md).
For every gene set, include the full original and used lists, species, missing
or changed genes, provenance and scoring method. Preserve software/function/
version, normalization and input layer, analysis unit and key parameters; separate
score calculation from display scaling and group averaging. State the z-score
axis/reference/aggregation or GSEA ranking statistic and contrast where applicable.

End with concrete reuse proposals when useful. A reading task does not authorize
editing installed skills. For requested skill changes, follow
[references/promotion.md](references/promotion.md); use available skill-authoring
tools or edit the target directly and run its relevant check. No external
skill-development plugin is required.

Use `paper-distill` in the companion biomed-skill package for general biomedical
reading, and `paper-evidence-map` for a claim–experiment–result graph.
