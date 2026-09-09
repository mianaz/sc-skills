---
name: omic-catalog
description: Use when finding, browsing, or downloading public omics datasets (GEO, SRA, CellxGene) into a local catalog. Hands accessions to sc-dataretrieval for loading.
---

# Omics Dataset Catalog

Backed by the `omic-catalog-mcp` server. Three modes — pick by user intent.

## When to use
- "find sc datasets for <topic>" / "what datasets exist for <paper/disease>" → DISCOVER
- "show me what's in my catalog" / "how many have raw counts" → BROWSE
- "download raw counts / processed data for <accession>" → DOWNLOAD

## When NOT to use
- Validating/loading an already-chosen accession into Seurat → use `sc-dataretrieval`
- Preprocessing CellRanger output → use `sc-preprocessing`

## Supported sources (Phase 1)
`supported`: GEO (GSE/GSM), SRA (SRP/PRJNA), CellxGene. Literature mining via PubMed + Europe PMC open-access full-text.
Everything else is flagged, never dropped: ArrayExpress / Single Cell Portal / HCA / Zenodo → `flagged_pending_adapter`; OMIX / GSA / EGA / dbGaP / lab sites / unknown URLs → `flagged_manual`.

## Mode: DISCOVER — database-first
1. `catalog_query(fts_query)` — check local catalog first (instant, offline).
2. If thin (<5 hits): `search_datasets(query, sources, organism, tissue, disease, year_from, year_to)` — hits GEO + CellxGene; auto-upserts.
3. `detect_duplicates(accession_list)` — present duplicate groups (DOI/alias/fuzzy-title) for the user to confirm; keep the best-metadata record, link the rest as aliases. **Never auto-merge.**
4. Present results table (see Output format).

## Mode: DISCOVER — literature-first
Topic-first: `search_literature(query, year_from, year_to)` → list papers.
Paper-first: `get_paper_fulltext(pmid_or_doi)`.
Then: read the returned `data_availability` + `methods` text → extract accessions/URLs → `catalog_upsert` each supported/pending entry, flag manual ones with the source URL + paper citation. Honour the no-scraping rule: only PMC/Europe PMC open-access full-text; abstract-only fallback otherwise.

## Mode: BROWSE
`catalog_query(fts_query, filters)` + `catalog_stats()` → formatted table with tier and status indicators.

## Mode: DOWNLOAD
1. `list_formats(accession)` → present tier options per the `sc-dataretrieval` tier table (1=raw h5, 2=filtered h5/mtx, 3=processed h5ad/rds, 4=fastq, 5=normalized-only).
2. Infer tier from phrasing ("raw counts" → 1/2; "processed" → 3) or ask.
3. `prepare_download(accession, format_tier, destination)`.
   - If `manual_note` is returned (OMIX/GSA/EGA/unknown): relay it; do not attempt download.
   - If `size_estimate_gb > 5`: a `slurm_script` is returned — offer to submit it via the `hpc-etiquette` skill rather than running interactively.
4. Hand the accession + file path + tier to the `sc-dataretrieval` skill for validation and Seurat loading.

## Output format (DISCOVER result table)
```
✓  GSE198638    human NASH liver    GEO          raw h5 + filtered    2023
✓  CXG-a3b2...  human NASH liver    CellxGene    h5ad                 2022
~  E-MTAB-9765  human NASH liver    ArrayExpress [adapter pending Phase 2]
!  OMIX002341   mouse NASH liver    OMIX         [manual: ngdc.cncb.ac.cn/omix/OMIX002341]
!  lab.github   NASH atlas          Unknown      [manual: extracted from doi:10.1016/...]
```
Legend: `✓` validated/supported · `~` pending adapter · `!` manual verification needed.

## Notes
- Set `NCBI_API_KEY` (and `NCBI_EMAIL`) in the environment for 10 req/s Entrez instead of 3.
- Duplicate fuzzy-title matches (rapidfuzz ≥0.85) can have false positives — always confirm with the user before merging.
- This skill stops at "here is the accession and its tier"; `sc-dataretrieval` owns everything after.
