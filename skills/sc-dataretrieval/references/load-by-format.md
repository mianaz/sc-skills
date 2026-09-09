# Loading deposited data into Seurat, by format

Each format maps to a pipeline entry point. **Only true RAW counts go through CellBender** (sc-preprocessing step 1). For everything else, load directly and set `seu$cellbender_applied <- FALSE` plus a methods.md note that ambient RNA was not removed. Keep raw counts in the `counts` slot/layer — downstream scVI and scDblFinder need them.

## Tier 1 — 10x mtx triplet (barcodes / features / matrix)
```r
# directory holds barcodes.tsv(.gz), features.tsv(.gz) (or genes.tsv), matrix.mtx(.gz)
counts <- Seurat::Read10X("sampleA/")             # gz handled automatically
seu <- Seurat::CreateSeuratObject(counts, project = "sampleA")
```
If files are prefixed (`GSM..._barcodes.tsv.gz`), rename to the bare 10x names in a per-sample dir first, or use ReadMtx:
```r
counts <- Seurat::ReadMtx(
  mtx = "GSM_matrix.mtx.gz",
  cells = "GSM_barcodes.tsv.gz",
  features = "GSM_features.tsv.gz")
```
- **raw/unfiltered** triplet → hand the dir/h5 to sc-preprocessing CellBender.
- **filtered** triplet → skip CellBender; flag ambient-RNA caveat.

## Tier 1 — 10x h5
```r
counts <- Seurat::Read10X_h5("sampleA.h5")        # raw_ or filtered_feature_bc_matrix.h5
seu <- Seurat::CreateSeuratObject(counts, project = "sampleA")
```
`raw_feature_bc_matrix.h5` → CellBender. `filtered_feature_bc_matrix.h5` → skip CellBender.

## Tier 2 — count table (csv / txt / tsv), genes × cells
```r
m <- as.matrix(data.table::fread("GSE_counts.csv.gz"), rownames = 1)  # genes in col 1
# verify orientation: rownames should be genes, colnames cells. Transpose if reversed:
# if (mean(grepl("^[ACGT]+", rownames(m))) > 0.5) m <- t(m)   # barcodes-as-rows -> transpose
seu <- Seurat::CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE), project = "GSE")
```
Decide raw-vs-filtered from cell count and the GEO docs; route accordingly. If values are non-integer it is normalized, not raw → tier 5.

## Tier 3 — processed object
Convert to Seurat and **prefer the authors' raw `counts` layer** if present (lets you rejoin at tier 1/2). See sc-conventions references/r-python-bridge.md for the converters.
```r
# h5ad -> Seurat (anndataR preferred; scConvert / SeuratDisk alternatives)
seu <- anndataR::read_h5ad("study.h5ad", as = "Seurat")
# if adata.raw / a counts layer exists, ensure it lands in the Seurat counts slot

# rds / h5Seurat
seu <- readRDS("study.rds")                       # may already be a Seurat object
seu <- SeuratDisk::LoadH5Seurat("study.h5Seurat")

# loom
seu <- SeuratDisk::Connect("study.loom", mode = "r") |> Seurat::as.Seurat()
```
Check what's inside before trusting it: `Seurat::Assays(seu)`, `slotNames`, `Layers(seu)` — is there a raw counts layer, or only normalized/scaled data? Record this.

## Tier 4 — reprocess from fastq
No usable counts deposited. Run CellRanger to produce `raw_feature_bc_matrix.h5`, then enter at tier 1 (→ CellBender). Match the authors' reference genome/transcriptome version (from provenance.md) where feasible.
```bash
cellranger count --id sampleA --transcriptome /refdata/refdata-gex-GRCh38-2024 \
  --fastqs fastqs/sampleA --sample sampleA --create-bam false
```

## Tier 5 — normalized matrix only (TPM/CPM/logged, no counts)
Load like tier 2, but **flag heavily**: CellBender, scDblFinder, scrublet, and scVI all expect counts and may be invalid. Set `seu$counts_are_normalized <- TRUE` and note in methods.md which downstream steps were skipped or adapted.

## After loading (all tiers)
- Set per-sample metadata (sample id, condition, accession) so a later merge in sc-integration can group correctly.
- Set the ambient/normalization flags above so the caveat travels with the object.
- Then proceed: raw → sc-preprocessing CellBender; filtered/processed → sc-preprocessing doublet flagging + QC (CellBender step skipped).
