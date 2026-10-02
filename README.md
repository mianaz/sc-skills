# sc-skills

Seurat-first single-cell analysis. Based on the live
[mianaz/sc-skills](https://github.com/mianaz/sc-skills/tree/0a980b6) repository,
retrieved 2026-10-02. This package keeps its workflow names and analysis content;
it moves general plotting and reproducibility to the companion `biomed-skill`.

## Workflows

| Stage | Skills |
|---|---|
| Entry and conventions | `single-cell`, `sc-conventions` |
| Data discovery and retrieval | `omic-catalog`, `sc-dataretrieval` |
| Core analysis | `sc-preprocessing`, `sc-integration`, `sc-annotation` |
| Sample-level condition tests | `sc-pseudobulk`, `sc-differential-abundance` |
| Downstream biology | `sc-trajectory`, `sc-cellchat`, `sc-grn`, `sc-spatial`, `sc-crispr`, `sc-target` |
| Single-cell literature methods | `sc-paper-distill` |

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
