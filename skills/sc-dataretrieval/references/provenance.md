# Capturing the authors' documented preprocessing

Public data comes with the authors' processing choices. Record them — they explain what state the matrix is in (raw vs filtered vs normalized), and decide what you should NOT redo. Append to methods.md alongside your own steps, clearly attributed to the original authors.

## Where the documentation lives
- **GEO**: the GSM sample page and GSE series page fields:
  - **Data processing** — alignment tool + reference genome/version, cell calling, filtering thresholds, normalization, batch correction, doublet handling, software versions. This is the primary source.
  - **Library strategy / Extraction protocol** — chemistry (10x 3'/5', version), library prep.
  - **Overall design / Summary** — sample/condition structure.
  - Pull these programmatically:
    ```r
    library(GEOquery)
    gsm <- getGEO("GSM1234567")
    Meta(gsm)$data_processing      # often a multi-line vector — keep all lines
    Meta(gsm)$extract_protocol_ch1
    Meta(gsm)$characteristics_ch1  # sample-level condition/tissue labels
    ```
- **ArrayExpress / BioStudies**: the **IDF** (protocols, overall design) and **SDRF** (per-sample → file mapping, processing) files; download and store them.
- **The paper + supplement**: Methods section usually has the authoritative pipeline; cite it.

## What to record in methods.md
- Accession(s) and per-sample file taken + its tier (from SKILL.md table).
- Alignment/quantification: tool (CellRanger/STARsolo/kallisto|bustools/alevin), reference genome + annotation version.
- Cell calling & filtering: how cells were called; nFeature/nCount/percent.mt cutoffs the authors applied.
- Normalization/transformation already applied (critical for tier 3/5 — tells you if `counts` are real counts).
- Doublet removal, ambient-RNA correction (e.g. CellBender/SoupX/DecontX), integration/batch method the authors used.
- Whether their object is raw, filtered, or normalized — and therefore which of YOUR steps to skip vs run.

## Reconciling their processing with this pipeline
- If authors already **filtered cells**: you cannot run CellBender (no empty droplets) — flag ambient RNA not removed; still flag doublets + QC.
- If authors already **removed doublets**: note it; running scDblFinder/scrublet again is optional and may double-flag.
- If authors already **normalized** (tier 5): treat as analysis-of-processed-data; do not pretend you have counts.
- If you intend to re-derive everything: get to raw counts (tier 1/2/4) and document that you deviate from the authors' pipeline and why.

Always make clear in methods.md which decisions are the **original authors'** vs **yours**.
