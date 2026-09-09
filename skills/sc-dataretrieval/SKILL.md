---
name: sc-dataretrieval
description: Use when loading a known public scRNA-seq accession (GEO, ArrayExpress, SRA) into Seurat after validating species/assay/tissue and choosing a deposited format tier.
---

# sc-dataretrieval (step 0 — get public data in)

## Overview
Turns a public accession (GSE…, E-MTAB…, SRP…) into an analysis-ready Seurat object **and** a record of how the authors processed it. The deposited format dictates where you enter the pipeline: only true RAW counts get CellBender; everything else skips ambient removal and must be flagged as such. Always also capture the authors' documented preprocessing into methods.md — see references/provenance.md.

**Validate before you trust the label.** A deposited series `title`/`gdsType` and any pre-existing curation field are hypotheses — confirm species, assay, tissue, and sample-vs-individual count from GSM-level metadata + the paper first (references/validation.md). Wrong-tissue, bulk-mislabeled-as-single-cell, and human-labeled-but-human+mouse accessions are common; there's no point format-tiering a series that turns out to be the wrong assay or organism.

## Format priority — take the highest tier available
Pick the best format the authors deposited; don't reprocess from fastq if usable counts already exist.

| Tier | Format | Pipeline entry | Notes |
|---|---|---|---|
| 1 | 10x **raw/unfiltered** counts — h5 (`raw_feature_bc_matrix.h5`) or mtx triplet | → CellBender → sc-preprocessing | only RAW gets ambient removal |
| 1 | 10x **filtered** counts — h5 (`filtered_feature_bc_matrix.h5`) or mtx triplet | → load → **skip CellBender** | flag: no ambient removal possible |
| 2 | raw count table — csv/txt/tsv (genes × cells) | → load → branch by raw-vs-filtered | infer filtered/raw from cell count & docs |
| 3 | processed object — h5ad / rds / h5Seurat / loom | → convert → carry in | preserves authors' normalization; prefer their raw `counts` layer if present |
| 4 | reprocess from **fastq** (SRA/ENA) | → CellRanger → tier-1 raw → CellBender | only when no usable counts exist; expensive |
| 5 | other normalized matrices (TPM/CPM/logged, no counts) | → load if usable | last resort; many count-based steps (CellBender, scVI, scDblFinder) won't apply — flag heavily |

**Always prefer counts over normalized values**, and within counts prefer raw over filtered. A tier-3 object that contains a raw `counts` layer can be treated as tier-1/2 once converted.

## Steps
1. **Validate the accession** → references/validation.md. Confirm species, assay/data_kind, tissue, and n_samples-vs-n_individuals from GSM-level metadata, and cross-check the publication. Record `validation_status` + a note. If it's the wrong tissue/assay/organism, stop here — don't download.
2. **Find & download** the data for the accession → references/sources.md (GEO, ArrayExpress/BioStudies, SRA/ENA).
3. **Identify the tier** of each deposited file using the table above; record what was available and what you chose.
4. **Load at the right entry point** → references/load-by-format.md. Raw 10x → hand to sc-preprocessing CellBender; filtered/processed → load directly and set a `cellbender_applied = FALSE` / ambient-RNA caveat flag.
5. **Capture provenance** — the authors' documented preprocessing (alignment ref, filtering, normalization, doublet handling) → references/provenance.md → methods.md.

Log the accession, its validation verdict, the per-sample file chosen + its tier, download date, checksums if given, and the documented preprocessing in methods.md.

## Relationship to omic-catalog
This skill takes **one known accession** and turns it into a loaded, provenance-tracked object. For *discovering* which datasets exist (searching/deduping/cataloguing across GEO/SRA/CellxGene into a local catalog), use **omic-catalog** first — it hands the chosen accession(s) here for validation and loading. If you already have the accession, start here.

## When NOT to use
Data already on disk as CellRanger output → straight to sc-preprocessing. Discovering/browsing candidate datasets before you have an accession → omic-catalog. Loading mechanics for a format you've already downloaded are in references/load-by-format.md; the rest of the pipeline is sc-preprocessing onward.
