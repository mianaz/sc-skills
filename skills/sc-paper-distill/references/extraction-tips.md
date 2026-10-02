# Extraction tips (finding the paper, reading code, mapping code→figures)

## Source order

Read supplied local PDFs, supplements and code first. Use supplied links next,
then retrieve missing public material. Ask only when missing information prevents
identifying the paper or the requested task.

## Locating the paper from a repo

When the user gives only a code URL:

1. Read the repo `README.md` — it usually names the paper, journal, and a DOI or a "cite this" block.
2. Check for a `CITATION.cff`, a DOI badge, or a `paper/` directory.
3. Look at the GitHub org/owner (e.g. `meta-cancer`) and repo description for the project name.
4. If still unresolved, search by title/first-author with the `paper-lookup` / `pubmed-research` / `literature-review` skills, or `gh search`.
5. Record the resolved DOI in the digest front-matter. If the paper is paywalled and you only have the abstract + figures, say so and work from captions.

## Full-text retrieval fallbacks

- Check PMID/PMCID first; many "paywalled" journals have OA mirrors in PMC.
- If publisher/preprint URLs are blocked (e.g. Cloudflare), use NCBI E-utilities:
	- `esearch`/`esummary` for PMID lookup
	- `efetch?db=pmc&id=<PMCID>&retmode=xml` for full-text XML
- Use Europe PMC `resultType=core` to discover `fullTextUrlList` and PMCID mappings.
- If only methods are accessible (no figures/repo), mark evidence level explicitly (`methods-only`).

## Pin the commit

Always record the exact commit you read: `gh api repos/<owner>/<repo>/commits/HEAD --jq '.sha'`. Snippets cited as `file:line` are only meaningful against a pinned SHA — the repo will drift.

## Reading the code efficiently

- Pull the file tree first (`gh api repos/<owner>/<repo>/git/trees/HEAD?recursive=1`) and read the *names* — numbered pipelines (`0-1`, `1a-2`, `2-2`) encode the analysis order.
- Read whole short scripts; for long ones, read the parts that set parameters, choose palettes/themes, or build the plotted object.
- Prefer `gh api .../contents/<path>` (or raw URLs) over cloning for a handful of files. For a deep dive, clone once into a scratch dir.
- For `.Rmd`/`.ipynb`, the plotting chunks are the figure code; the prose chunks often state the figure's intent.

## Mapping code → figures/methods

The hard part: papers rarely label code by figure number. Infer the mapping from:
- File/chunk names (`subclus_marker_plot` → the marker dotplot/heatmap panels).
- Object names that match caption nouns ("Hepatocyte", "TAM", "malignant").
- The pipeline position (preprocessing scripts → QC supplementary panels; subclustering → the per-compartment UMAPs).
State the mapping as your inference, with confidence, rather than asserting it as fact.

## Handling partial or missing figure code (the common case)

Most repos ship analysis code and only some plotting code. When figure code is absent:
- Still write a **reconstructed template** in §3 that reproduces the *style* implied by the caption, using house tools (scop / tidyplots / ggplot2 / ComplexHeatmap per `scientific-plotting`).
- Mark it clearly: `Code: caption-only — no code provided; template reconstructed in house style`.
- Do not invent parameter values the paper didn't state; choose sensible house defaults and flag them as choices.

## Marker/signature extraction discipline

- Extract markers with full context, not as bare gene lists:
	- experiment (flow/scRNA/CITE/Visium/IF/qPCR/WB), species, and biological context
	- what claim they support (identity, zonation, perturbation response, niche signaling)
- Store each marker/signature-context unit as one row in the output `omic_marker_db.csv`.
- If a signature is referenced but not fully enumerated in text, log it with a note like `list not fully enumerated in manuscript text`.

## What counts as "notable" (don't catalog everything)

Catalog a method/figure only if it is **reusable or instructive**: a new or non-obvious method, a parameter choice worth copying, or a figure whose style we'd want to reproduce. Skip boilerplate (standard `NormalizeData`/`RunUMAP` with default params) unless the paper does something deliberate with it. A tight digest of 4–8 figures and 3–6 methods beats an exhaustive dump.
