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

Curation rules:
- One row per marker/signature-context pair.
- Do not infer markers that are not explicitly present in paper text/figures/tables.
- If a signature is named but not fully listed, store the named signature and add a caveat in `notes`.
- Keep rows factual; interpretation belongs in digest files.
