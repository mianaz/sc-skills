# Integration QC

## Default — visual
UMAP colored by batch/sample vs. by cell type / canonical markers. Good mixing of
batches WITHOUT collapsing biological structure indicates acceptable integration.
Render via scientific-plotting (CellDimPlot), dual-version + source data.

## Optional — quantitative (when rigor needed)
scIB metrics (batch: iLISI, kBET; bio: cLISI, ARI, NMI). Document the score table
in methods.md. Only run when reviewers/decisions require it.
