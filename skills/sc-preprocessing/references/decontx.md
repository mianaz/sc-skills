# decontX — R-native ambient-RNA removal (alternative to CellBender)

Use when you do NOT have CellRanger RAW output (only the filtered matrix), or want a
pure-R per-sample pipeline. decontX estimates and removes ambient contamination from a
filtered count matrix; CellBender needs the raw empty-droplet tail, decontX does not.

```r
library(SingleCellExperiment)
library(decontX)

# dat: genes x cells filtered counts (e.g. from Read10X_h5(...)[["Gene Expression"]])
sce <- SingleCellExperiment(list(counts = dat))

# drop truly empty barcodes first
lib.sizes <- colSums(counts(sce))
sce <- sce[, lib.sizes > 0]

sce <- decontX(sce)                       # adds decontXcounts(sce) + sce$decontX_contamination
seu <- Seurat::CreateSeuratObject(
  counts = round(decontXcounts(sce)),     # round the corrected counts back to integers
  project = sample_id, min.cells = 1, min.features = 1)
```

Run **per sample, before merging**. Log mean `decontX_contamination` per sample in methods.md.
Doublet flagging (scDblFinder) can run on the same `sce` in the same loop — see doublets.md.

## Decision: CellBender vs decontX
- True RAW matrix on disk, GPU available → CellBender (cellbender.md).
- Only filtered matrix, or R-only environment → decontX (here).
- Neither possible → skip correction, flag ambient, note it in methods.md.
