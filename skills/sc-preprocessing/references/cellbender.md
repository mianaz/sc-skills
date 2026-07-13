# CellBender (per-sample, on CellRanger RAW output)

Run on `raw_feature_bc_matrix.h5` (the RAW matrix — not filtered), once per library.
```bash
cellbender remove-background \
  --input  raw_feature_bc_matrix.h5 \
  --output sampleA_cellbender.h5 \
  --cuda \
  --expected-cells <N> \
  --total-droplets-included <M>   # > expected-cells; include the empty-droplet tail
  # --epochs left at default unless convergence/QC indicates otherwise
```
Parameters are dataset-specific: set `--expected-cells` from the CellRanger estimate and `--total-droplets-included` to capture the ambient tail; record chosen values + rationale in methods.md. Output `sampleA_cellbender_filtered.h5` holds denoised counts.

Load into Seurat:
```r
counts <- Seurat::Read10X_h5("sampleA_cellbender_filtered.h5")
seu <- Seurat::CreateSeuratObject(counts, project = "sampleA")
```
