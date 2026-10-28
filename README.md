# sc-skills

Seurat-first single-cell analysis, from public data retrieval through annotation,
sample-level tests and downstream biology. General plotting and reproducibility
are provided by the companion [biomed-skill](https://github.com/mianaz/biomed-skills).

## Workflows

| Stage | Skills |
|---|---|
| Entry and conventions | `single-cell`, `sc-conventions` |
| Data discovery and retrieval | `omic-catalog`, `sc-dataretrieval` |
| Core analysis | `sc-preprocessing`, `sc-integration`, `sc-annotation` |
| Sample-level condition tests | `sc-pseudobulk`, `sc-differential-abundance` |
| Downstream biology | `sc-trajectory`, `sc-cellchat`, `sc-grn`, `sc-spatial`, `sc-crispr`, `sc-target` |
| Single-cell literature methods | `sc-paper-distill` |

## Gene sets

Prefer UCell for per-cell gene-set scoring. Use cited published signatures/data
or unbiased DE signatures from the user's data. Never invent a gene set heuristically.
See the [scoring and reuse rules](skills/sc-conventions/references/gene-sets.md)
and the [curated resource guide](skills/sc-conventions/references/gene-set-resources.md).

Keep reusable literature sets in `SC_GENE_SET_HOME`, defaulting to
`~/.local/share/sc-skills/gene-sets/`. The suite includes a source-checked
[12-chemokine human signature](skills/sc-conventions/assets/gene-sets/messina_2012_tls_12chemokine_human_v1.json)
as a starter record. Third-party records retain the source terms recorded with each set.

## Shared dependency

Install [`biomed-skill`](https://github.com/mianaz/biomed-skills) alongside this package for `scientific-plotting` and
`scientific-reproducibility`. Stage skills refer to those shared names; they have
one implementation in biomed-skill. `sc-paper-distill` stays specific to omics
method curation. General paper reading and evidence maps belong in biomed-skill.

For Claude Code:

```bash
claude plugin marketplace add mianaz/biomed-skills
claude plugin install biomed-skill@biomed-skill
claude plugin marketplace add mianaz/sc-skills
claude plugin install sc-skills@sc-skills
```

For Codex or Cursor, copy both packages' skill folders into the same agent skill
directory, preserving existing local edits before replacement.

Example: “Use $single-cell to analyze these samples, then use
$scientific-plotting for the requested figures.”

MIT; upstream attribution is retained in LICENSE.
