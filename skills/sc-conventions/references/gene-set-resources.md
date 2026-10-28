# Gene-set and cell-marker resources

Curated on 2026-10-28. Use the [selection and storage rules](gene-sets.md).
For a ready-to-score set, start with a complete published signature or MSigDB
collection. For cell identity, use a marker database with a species/tissue query.
Choose state/pathway signatures for state questions and cell-type signatures for
identity questions.

| Resource | Best use | Retrieval and citation |
|---|---|---|
| [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) | Published pathways and programs; human C8 and mouse M8 for cell-type signatures; H/MH for hallmarks; C7/M7 for immune states | Download a named set or collection in GMT, or use `msigdbr`. Record database species, release, collection, set ID and original study. Follow the [citation guidance](https://www.gsea-msigdb.org/gsea/msigdb/) and cite the set's source paper. |
| [CellMarker 2.0 / 3.0](https://bio-bigdata.hrbmu.edu.cn/CellMarker/) | Human/mouse cell identity with tissue, disease and publication context | Use the site's Download page or a saved query. Cite [Hu et al., 2023, 2.0](https://doi.org/10.1093/nar/gkac947) or [Hu et al., 2026, 3.0](https://doi.org/10.1093/nar/gkag1040), matching the release used, plus the source studies. |
| [PanglaoDB](https://panglaodb.se/markers.html) | Human/mouse marker lookup and comparison across cell types | The marker page provides an export; the [FAQ](https://panglaodb.se/faq.html) describes retrieval. Retain species tags, cell type, tissue and available specificity/evidence fields. Cite [Franzén et al., 2019](https://doi.org/10.1093/database/baz046). |
| [Cell Taxonomy](https://ngdc.cncb.ac.cn/celltaxonomy/) | Cell types across species, tissues and conditions, with supporting studies | Download `Cell_Taxonomy_resource.txt` from the [official download page](https://ngdc.cncb.ac.cn/celltaxonomy/download). Preserve its `Source` field and publication records. Cite [Jiang et al., 2023](https://doi.org/10.1093/nar/gkac816). |
| [Human Protein Atlas: Single Cell](https://www.proteinatlas.org/humanproteome/single+cell) | Human cell-type expression and specificity, with tissue/protein context | Use the [download page](https://www.proteinatlas.org/about/download) or export a documented specificity query. Retain release, tissue, cell type and RNA/protein evidence type. Cite [Karlsson et al., 2021](https://doi.org/10.1126/sciadv.abh2169) and the versioned data page per [HPA guidance](https://www.proteinatlas.org/about/licence). |

## Resource-specific choices

- **MSigDB:** the site lists human and mouse release 2026.2 at this curation date.
  Pin the release actually retrieved, including `db_version` when using `msigdbr`.
  Prefer native-species collections where suitable; retain any ortholog mapping.
  Web downloads require registration. Preserve the source terms for each collection.
- **CellMarker:** 3.0 was published on 2026-10-27 and the main site now identifies
  as 3.0. It adds computational predictions alongside literature-curated records.
  Prefer the literature-curated evidence for established markers. Save the evidence
  class, method and dataset for predicted records. Preserve older 2.0 files as
  versioned sources; a fresh download from the current site needs its own version.
- **PanglaoDB:** save the export date or filename. Its human/mouse tags need to be
  resolved before using the markers in a species-specific score. Inspect each
  record's source support and cell-type context.
- **Cell Taxonomy:** its [help page](https://ngdc.cncb.ac.cn/celltaxonomy/help)
  documents integration of CellMarker, PanglaoDB and other sources. Deduplicate
  support by original study and marker-context record when combining searches.
- **HPA:** select genes by a documented specificity classification or quantitative
  rule applied to the exported table. Keep the exact query and result. Its
  normal-tissue reference is useful for lineage identity; use disease-matched
  studies when defining a disease-state program.

## Add a resource or literature set

Add a resource here when it has a clear use, a primary citation and a retrievable
table, API or author repository. Record the retrieval route and relevant evidence
fields. Cache only the release needed for the task in the shared library.

For a paper-derived signature, retrieve the complete membership from the source,
save the table or extraction record, and add a versioned JSON definition to
`SC_GENE_SET_HOME/sets/` using the [shared schema](gene-sets.md). Reuse it on later
projects. Keep gene-level citations when a database selection combines studies.
