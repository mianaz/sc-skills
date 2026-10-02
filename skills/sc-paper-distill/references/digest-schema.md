# Digest schema (fixed template — fill one per paper)

Save as `papers/<year>-<shortname>.md` (e.g. `papers/2026-scPLC.md`). Use this exact section order so every digest is comparable. Omit a section only if it genuinely does not apply, and say so explicitly rather than deleting the heading.

````markdown
---
paper: "<short name> — <full title>"
authors: "<first author> et al."
doi: "<doi>"
journal: "<journal>"
year: <year>
repo: "<repo url> @ <commit sha>"
data_accession: "<GEO/SRA/CELLxGENE/EGA ids, if any>"
distilled: <YYYY-MM-DD>
tags: [<organism>, <tissue/disease>, <modality>, <method keywords>]
---

# <short name> — distillation

## 1. Overview
- **What it is & why notable:** 2–4 sentences. Why did we pick this paper — pleasing figures, new method, or both?
- **Data:** organism, tissue/disease, platform (10x 3'/5', Visium, …), n samples / n cells, modalities.
- **Code availability:** what the repo actually ships (analysis only? plotting? both?) vs. what is caption-only / absent. Be specific by file.

## 2. Methods catalog
One block per notable method. Number them M1, M2, …

### M1. <method / analysis name>
- **Purpose:** why used in this paper.
- **How used:** key parameters and choices (the decisions that matter for reuse).
- **On what data:** which step / which cells / which comparison.
- **Tools / packages:** with versions if stated.
- **Code:** `path/in/repo.R:LL-LL` — short snippet or template (or "not provided").
- **Interpretation:** what the result showed / the biological or technical claim.
- **Reusable?** yes/no → target skill (e.g. sc-grn, sc-preprocessing) + one-line what to add.

## 3. Figure catalog
One block per notable figure or panel. Number them F1, F2, …

### F1. Fig <n><panel> — <one-line description>
- **Type:** UMAP / dotplot / heatmap / stacked-bar composition / volcano / river / spatial / …
- **What it communicates:** the message the panel carries.
- **Data shown:** what's on each axis / what the color & size encode.
- **Visual style:** palette (named if identifiable), key aesthetic choices (point size, ordering, white borders, faceting, label strategy), layout.
- **Code:** repo file `path:LL` OR "caption-only — no code provided".
- **Reusable template:** the minimal snippet to reproduce the *style* on the target data (write one even if the paper gave no code — mark it as reconstructed).
- **Interpretation:** what the reader is meant to conclude.
- **Reusable?** yes/no → target skill (scientific-plotting / sc-conventions) + one-line what to add.

## 4. Style observations
Cross-cutting aesthetic patterns seen across this paper's figures: palette family, fonts, sizing/aspect conventions, recurring layout idioms, how they annotate stats. These feed sc-conventions, not individual recipes.

## 5. Promotion Proposals
The actionable output. One row per candidate edit to the sc suite. The user approves/edits/rejects each.

| # | target skill | kind | summary | source |
|---|---|---|---|---|
| P1 | scientific-plotting | recipe | <what recipe/template to add> | F1 |
| P2 | sc-conventions | requirement | <what rule/pattern to add or strengthen> | §4 |
| P3 | sc-grn | method | <what method recipe to add> | M2 |

`kind` ∈ { recipe, requirement, method, palette, fix }.

## 6. Directed distill tables (when requested)

### 6.1 Paper info table

Include one markdown table row with:

| title | authors | journal | year | doi | main findings summary |
|---|---|---|---:|---|---|
| ... | ... | ... | ... | ... | ... |

And persist the same row in `references/paper_catalog.csv`.

### 6.2 Marker/signature curation table

Include a marker/signature table with at least these columns:

| marker_or_signature | marker_type | experiment | species | context | claim_use | evidence_source | citation |
|---|---|---|---|---|---|---|---|
| ... | ... | ... | ... | ... | ... | ... | ... |

Then append rows to `references/omic_marker_db.csv`.

Use one row per marker/signature-context pair (do not merge unrelated contexts into one row).

## 7. Evidence level tags

For each method/figure/code claim, add one evidence-level note when coverage is partial:

- `full-text`
- `methods-only`
- `captions-only`
- `abstract-only`
- `repo-only`
````

## Filling rules

- **Number everything** (M1/F1/P1) so §5 can cite exact sources and approvals are unambiguous.
- **Snippets are templates, not transcriptions.** Strip paper-specific paths/object names down to the reusable shape, but keep the parameter choices that define the style/behavior.
- **A figure with no code still gets a reconstructed template** in §3 — that is often the most valuable output, since it converts a caption into runnable house-style code.
- **§5 must be self-contained:** a reader should be able to act on each row without re-reading §2/§3, but the `source` column lets them.
