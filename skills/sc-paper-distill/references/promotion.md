# Promotion (digest §5 → the sc suite)

The digest is paper-specific and append-only. The suite is cross-paper and curated. Promotion is the controlled bridge between them.

## The gate

After writing the digest, present the §5 Promotion Proposals table. The user approves / edits / rejects each row. **Do not edit any existing sc skill before approval.** Rejected rows stay in the digest (still useful as a record) but never touch the suite.

## Routing table

| `kind` | Target | What it becomes |
|---|---|---|
| `requirement` (universal) | `scientific-reproducibility` | A publication-grade output/reproducibility rule that is NOT single-cell-specific (figure contract, logging, dpi, stats). |
| `requirement` (single-cell) | `sc-conventions` (a new `patterns.md` if no existing home) | An sc-specific rule or recurring style pattern (Seurat handling, cell-type ordering, atlas conventions). |
| `palette` | `sc-conventions/references/palettes.md` | A named palette / color convention, with when-to-use. |
| `recipe` | `scientific-plotting/references/*` (the R plotting skill; sc-object recipes stay scop-flavored, general figure recipes are not sc-specific) | A reusable figure template (new reference file, or appended to the matching one — volcano, heatmaps, scop-plots, …). |
| `method` | the matching method skill — `sc-*` for single-cell, or a non-sc skill (`bulk-rnaseq`, `pydeseq2`, …) when the recipe belongs there | A reusable analysis recipe at the right pipeline stage. |
| `fix` | wherever the error lives | A correction to existing skill content the paper revealed. |

When a recipe doesn't fit any existing `scientific-plotting` reference, create a new `references/<topic>.md` and add a row to the SKILL.md tool router — don't bloat an unrelated file.

## Apply promotions carefully

Preserve the target skill's quality bar when applying a promotion (baseline → small edit → sanity-check). If you have a skill-authoring helper installed (e.g. `superpowers:writing-skills`), route the edit through it; otherwise edit the target `SKILL.md` directly. Batch the approved rows for one paper into a single pass when they touch the same skill.

## Provenance back-link (mandatory)

Every promoted item carries a back-link to its source digest so the suite stays traceable:

```r
# Stacked composition with significance brackets.
# (source: sc-paper-distill/papers/2024-example.md F4)
```

or in prose references:

> Use white tile borders on square heatmaps (source: papers/2024-example.md §4).

When a later paper strengthens an existing rule, append its source rather than replacing the first:
`(source: papers/2024-example.md §4; papers/2025-example2.md §4)`.

## Dedup before adding

Before adding a `requirement` or `palette`, check whether `sc-conventions` already states it. If so, strengthen/extend the existing line (and append the new source) instead of duplicating. The corpus is where every instance is recorded; the suite holds each rule once.

## Marker/signature promotions

Marker/signature rows in `references/omic_marker_db.csv` are evidence inventory, not automatic skill rules. Promote only when:

- the same marker/signature pattern appears across multiple papers, or
- a marker panel is central and well-validated for a narrow domain (for example, liver KC identity).

When promoted, route to:

- `sc-annotation` for cell-identity marker panels
- `sc-spatial` for zonation/niche signatures
- `sc-target` for perturbation/axis-linked signatures

Always keep the promoted snippet linked back to paper IDs in the marker DB and digest sources.
