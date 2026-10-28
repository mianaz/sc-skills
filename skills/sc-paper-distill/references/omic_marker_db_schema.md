# omic_marker_db.csv schema

`omic_marker_db.csv` is a marker/signature evidence table created in the user's
research output directory. It is not bundled with the skill.

Columns:
- `paper_id`: digest identifier, e.g. `2022-hepatic-macrophage-niches`
- `marker_or_signature`: gene/protein marker name or signature label
- `marker_type`: `gene/protein`, `signature`, `ligand signature`, `modeling marker set`, etc.
- `experiment`: assay where used (flow, scRNA, CITE-seq, Visium, IF, qPCR, WB, etc.)
- `species`: species scope (`mouse`, `human`, `mouse|human`, ...)
- `context`: biological/technical context (cell identity, zonation, niche sender set, perturbation)
- `claim_use`: what claim this marker supports
- `evidence_source`: where to verify (figure, STAR Methods, table)
- `citation`: full citation anchor (DOI at minimum)
- `notes`: caveats and details
- `added_on`: YYYY-MM-DD
- `gene_set_id`: shared-library set ID for a signature; empty for individual markers
- `gene_set_record`: path to its versioned JSON definition in the [shared library](../../sc-conventions/references/gene-sets.md)
- `membership_status`: `complete` or `incomplete` for a signature

Curation rules:
- One row per marker/signature-context pair.
- Do not infer markers that are not explicitly present in paper text/figures/tables.
- If a signature is named but not fully listed, save an `incomplete` library record
  with the source locator and the missing information in `notes`; exclude it from scoring.
- Store complete gene membership in the shared-library record and link it here.
  Reuse the same record across digests rather than copying partial gene lists.
- Keep rows factual; interpretation belongs in digest files.
