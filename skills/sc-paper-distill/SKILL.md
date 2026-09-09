---
name: sc-paper-distill
description: Distill a single-cell / spatial / omics paper plus its code into an append-only digest and promotion proposals for the sc skill suite.
disable-model-invocation: true
---
# Paper Distillation (paper + code → reusable skill knowledge)

## Purpose

This is the **producer** half of a two-part system for growing the sc skill suite from the literature.

- **Producer (this skill):** turn one paper + its code into an append-only, structured **digest** of paper-specific findings — what method, for what purpose, how used, on what data, what visualization, the interpretation — plus extracted reusable code snippets and provenance. Digests are never overwritten; the corpus is the lossless record.
- **Consumer (the sc suite):** the digest ends with **Promotion Proposals**. Approved proposals are applied via `superpowers:writing-skills` into `sc-conventions` (publication-grade requirements / style patterns), `scientific-plotting` (reusable figure recipes), and the method skills (`sc-grn`, `sc-annotation`, `sc-preprocessing`, …). The suite stays lean, curated, and deduplicated.

Keep these halves separate. The digest is additive and paper-specific; the suite is curated and cross-paper.

## Workflow

Read `references/digest-schema.md` and `references/extraction-tips.md` before starting. Do not work from memory of a paper — read the actual text and the actual code.

1. **Asset check (ask first).** Before searching the web, explicitly ask whether the user already has local readable assets:
	- full text PDF
	- supplementary data/tables
	- analysis code repo or local code folder
	If yes, read local assets first.
	If not, ask whether the user has ready links.
	If not, then search and resolve links yourself.
2. **Gather inputs.** You need the paper (URL / DOI / PDF) and the code repo URL. If only a repo is given, locate the paper first (README, DOI badge, repo → publication). See `references/extraction-tips.md`.
3. **Read the paper.** Relevant methods sections, the captions of the notable figures, and the results/interpretation tied to those figures. Capture *why* each method/figure exists, not just *what* it is.
4. **Read the code.** Clone/browse the repo at a pinned commit. Map code files → figures/methods. Extract reusable snippets and templates. Explicitly record what is **missing** (papers often ship analysis code but not figure code, or only some panels).
5. **Write the digest.** Fill the fixed schema → `papers/<year>-<shortname>.md`. One block per notable method (§2) and per notable figure (§3); roll up reusable items into the §5 Promotion Proposals table.
6. **Directed distill add-ons (required when requested).**
	- Add/update a row in `references/paper_catalog.csv` with title, authors, journal, year, DOI, accession, code status, and one-line findings summary.
	- Add a marker/signature curation block in the digest and append marker rows to `references/omic_marker_db.csv`.
	- Marker rows must include marker/signature, experiment, species, context, claim-use, citation, and provenance source.
	- If full text or code is unavailable, state `not available` explicitly and tag evidence level (abstract-only, captions-only, methods-only, full-text).
7. **Present Promotion Proposals.** Show the §5 table. The user approves / edits / rejects each row. This is a hard gate — do not edit the suite before approval.
8. **Apply approved promotions.** Invoke `superpowers:writing-skills` and route each approved row to its target skill, each carrying a provenance back-link `(source: sc-paper-distill/papers/<year>-<shortname>.md F3)`. See `references/promotion.md`.

## References

| Need | Reference |
|---|---|
| The fixed digest template every paper is analyzed against | `references/digest-schema.md` |
| Finding the paper from a repo, reading partial/missing code, mapping code→figures | `references/extraction-tips.md` |
| Routing approved proposals to the right skill + the provenance back-link convention | `references/promotion.md` |
| Persistent paper-level bibliography table | `references/paper_catalog.csv` |
| Persistent omics marker/signature database | `references/omic_marker_db.csv` |
| Marker/signature DB schema and curation rules | `references/omic_marker_db_schema.md` |
| Accumulated per-paper digests (the corpus) | `papers/*.md` |

## Non-Negotiables

- **One fixed schema, every paper.** Consistency is the point — it makes the corpus comparable and promotions mechanical.
- **Provenance always.** Every snippet records repo URL + commit SHA + `file:line`; every digest records the paper DOI/citation.
- **Note what's missing.** "Figure code not provided — reconstructed from caption" is a finding, not a gap to hide.
- **Append-only corpus.** Never overwrite a prior digest; add a new one (or a dated revision section) instead.
- **Approval gate before the suite changes.** §5 proposals are candidates until the user approves them.
- **Promotions go through `superpowers:writing-skills`** so existing sc skills keep their TDD/quality bar.
