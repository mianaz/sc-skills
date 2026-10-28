# Gene sets: selection, scoring and reuse

## Select an evidence-based set

**Prefer UCell over Seurat `AddModuleScore` for per-cell gene-set scoring.** Use
`AddModuleScore` when the user requests it or when reproducing an analysis that
requires it; record that choice. Retain a method's native scoring model when its
weights or statistics define the analysis, such as SCENIC/AUCell or scDRS.

Before scoring, search the shared library below, then the
[curated resources](gene-set-resources.md). Choose a set that matches the species,
tissue, cell population and biological question. Use one of these sources:

- A complete published signature from a paper, supplement, author code or
  versioned database. Cite the original study and the exact table, sheet, column,
  set ID or code location. For a published dataset, retain its accession and the
  reproducible extraction or DE rule as well as the citation.
- A signature derived from unbiased DE in the user's own data. Define the
  contrast, tested genes, model, FDR/effect-size filters, direction and any rank or
  size limit before examining the score pattern. Apply the rule to the full DE
  table and retain that table and the derivation code. Use replicate-aware DE for
  condition contrasts; label cluster-marker DE as such.

**Never make up a gene set heuristically.** Do not turn familiar genes, gene-name
patterns, narrative examples or a few highlighted paper genes into a signature.
Use the full source list or a recorded, systematic selection from a source table.
If membership cannot be recovered, save an `incomplete` record and find another
source or derive a DE signature. Incomplete records are excluded from scoring.

For marker databases, save the query and all records selected by species, tissue,
cell type, condition and evidence class. Preserve the original supporting studies.
Keep independently derived signatures separate; combining them requires an explicit
rule and a new derived-set record. Repeated citations across databases share provenance.

## Preserve the definition

- Keep the original gene IDs, order, signs and any weights. Save mapped symbols,
  duplicate handling, genes present in the input and missing genes separately.
- For another species, use a versioned orthology resource and save the mapping.
  Record unmapped and one-to-many genes. Never replace a missing ortholog with a
  related gene or infer orthology by changing letter case.
- Store the complete original set independently of each dataset's measured genes.
  Use a fixed definition and missing-gene policy across datasets being compared.
  Report coverage; assess whether the measured genes still cover the intended program.
- A short marker display panel retains its role as an annotation aid. A scored
  signature needs its own complete source membership.

Score mechanics, input layers, signed sets and output records are in the
[UCell recipe](../../sc-annotation/references/module-scores.md). For a signature
derived from the scored data, describe the result as discovery; use held-out
samples or an external cohort when testing its generalization.

## One shared library

Use `SC_GENE_SET_HOME` if configured; otherwise use
`~/.local/share/sc-skills/gene-sets/`. This user-owned library persists across
projects and plugin updates. Search it before downloading or extracting a set.
Also check the suite's [bundled published sets](../assets/gene-sets/).

```text
gene-sets/
  sets/<set_id>.json
  sources/<source_id>/<original table, supplement or extraction code>
```

Save useful literature sets here as they are extracted. Reuse a matching record
by source, source set ID, species and version. Add a new version for changed
membership or derivation; retain the earlier record. Keep private DE tables and
unpublished signatures in this user library or their research project.

Each JSON record contains:

| Fields | Content |
|---|---|
| `set_id`, `name`, `version`, `status` | Stable name/version; `complete` or `incomplete` |
| `source_type`, `species`, `gene_id_type`, `context` | Published signature, database selection or DE derivation; biological scope |
| `citation`, `source_url`, `source_locator`, `source_version`, `retrieved_on` | DOI/PMID or dataset accession, exact membership location, release/date |
| `genes_up`, `genes_down`, `weights` | Complete source membership and direction; empty arrays for unused directions; weights when present |
| `original_scoring`, `selection_rule` | Published method and exact extraction/filter rule |
| `source_files`, `source_license`, `notes` | Relative source paths when saved, source terms, relevant context |

For database selections, add `supporting_citations` and the saved query. For DE
sets, add `derivation` with the dataset, contrast, design, method, tested universe,
filter/ranking rule, discovery samples, full-results path and script path. For
mapped sets, add `parent_set_id` and the orthology resource/version/mapping path.

The [12-chemokine human set](../assets/gene-sets/messina_2012_tls_12chemokine_human_v1.json)
is a complete literature seed with a source locator and original scoring method.
Copy a bundled record into the shared library when first used. Save third-party
source files under their own terms; the suite's MIT license covers its instructions.

For each analysis, copy the selected JSON record to `source_data/gene_sets/` and
save the original-to-used gene table, coverage, score parameters and scores there.
This project snapshot fixes the version used while the shared library grows.

## Earlier literature leads

The earlier suite cited [TLS markers](https://doi.org/10.1126/science.adz2742) and
[Courau TME programs](https://doi.org/10.1038/s41590-026-02505-7). Recover their
complete source tables before adding them to the library; the former short code
lists are not retained as scoring definitions.
