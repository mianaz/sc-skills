# Terminology ledger — one name for one thing

The naming counterpart to the cell-type *order/color* rule in `palettes.md`. That rule
keeps a category in the same order and color across panels; this one keeps it under the
same **name** across every script, object column, figure, caption, and `methods.md`
entry. A label that drifts ("Kupffer cell" → "KC" → "Mono-Mac", `CD8 T` → `CD8+ T` →
`T_CD8`) forces the reader to re-learn it and reads as careless work.
(adapted from Yuan1z0825/nature-skills `_shared/core/terminology-ledger.md`.)

## 1. Build the ledger on first contact

When you first assign labels or wire up an analysis, extract every recurring term into a
canonical ledger **before** producing figures or tables:

- cell-type labels — coarse compartments and fine subtypes (and the cluster IDs they map to)
- gene / protein symbols — canonical HGNC (human) / MGI (mouse) casing; run through
  `HGNChelper` (see sc-preprocessing) rather than hand-casing
- sample / donor / condition / batch / timepoint names
- metric and assay names, units, statistical symbols
- abbreviations, each with its one first-use expansion

## 2. Store it as the single source of truth

Keep the ledger where the pipeline already reads from, not in prose:

- the canonical cell-type **`levels` vector** already used for order/palette
  (`ct_lv` in `palettes.md`) is the cell-type half of the ledger — reuse it, don't fork it
- a small `labels.csv` (cluster → coarse → fine → canonical label) under the project
- a terms block in `methods.md` (see scientific-reproducibility) for samples/metrics/abbreviations

## 3. Lock and enforce

- Use only canonical forms in every output — object columns, legends, dotplot axes,
  captions, source-data headers, filenames. Do not introduce synonyms to vary wording.
- Define each abbreviation once, then use the short form.
- If the user renames a term later, change **every** occurrence (columns, saved objects,
  figures, source data, methods.md) and update the ledger — not just the current script.
- Flag collisions explicitly: the same population under two names, or one name reused for
  two different populations.

## 4. Do not invent names

Do not coin a cell-type label (or rename a gene/sample) that the evidence does not
support. If a cluster's identity is unresolved, keep the cluster ID or a marker-based
provisional tag and flag it — never fill the gap with a guessed name. This is the naming
side of `sc-annotation`'s marker-verification rule: a label is a claim, and every claim
needs canonical-marker support with a citation.
