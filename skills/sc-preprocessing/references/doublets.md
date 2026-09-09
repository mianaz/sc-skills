# Doublet detection — FLAG, do not remove

Run per-sample, before merging.

## scDblFinder (R)
```r
library(scDblFinder)
sce <- scDblFinder(Seurat::as.SingleCellExperiment(seu))
seu$scDblFinder.score <- sce$scDblFinder.score
seu$scDblFinder.class <- sce$scDblFinder.class            # "singlet"/"doublet"
```

## scrublet (Python) — carry scores back
```python
import scrublet as scr
import pandas as pd
scrub = scr.Scrublet(counts_matrix)         # cells x genes (raw counts)
scores, predicted = scrub.scrub_doublets()  # predicted is boolean
pd.DataFrame(
    {"score": scores,
     "call": ["doublet" if p else "singlet" for p in predicted]},
    index=barcodes,                          # cell barcodes, aligned to counts_matrix rows
).to_csv("sampleA_scrublet.csv")             # -> columns: <barcode>,score,call
```
```r
sc <- read.csv("sampleA_scrublet.csv", row.names = 1)     # barcode,score,call
# align by barcode; mismatched barcodes (e.g. trailing "-1") yield NA -> fix before use
seu$scrublet.score <- sc[colnames(seu), "score"]
seu$scrublet.call  <- sc[colnames(seu), "call"]
```

## Consensus columns (nothing removed)
```r
seu$doublet_any <- (seu$scDblFinder.class == "doublet") |
                   (seu$scrublet.call == "doublet")
seu$doublet_consensus <- (seu$scDblFinder.class == "doublet") &
                         (seu$scrublet.call == "doublet")
```

## Bioconductor alternative (scran / DropletUtils)

A self-contained alternative stack for groups already in the SCE/Bioconductor world. Not the
house default (above), but a vetted option. (source: sc-paper-distill/papers/2022-scPLC.md M1, M3 — Xue et al. Nature 2022)

```r
# Empty droplets: data-driven lower bound (2nd-smallest total UMI) instead of a fixed number.
b <- sort(colSums(DropletUtils::... <- SingleCellExperiment::counts(sce)))[2]
e.out <- DropletUtils::emptyDrops(SingleCellExperiment::counts(sce), lower = b)
sce   <- sce[, which(e.out$FDR <= 0.01)]                 # keep FDR <= 0.01

# Doublets at the CLUSTER level: drop clusters whose median score exceeds median + 3*MAD.
scores  <- scran::doubletCells(sce)
lscore  <- log10(scores + 1)
cutoff  <- median(lscore) + 3 * mad(lscore)
# cluster (e.g. louvain on a KNN graph), then:
drop    <- tapply(lscore, clusters, median) > cutoff
sce     <- sce[, !clusters %in% names(drop)[drop]]
```

Trade-off vs. the house route: this *removes* doublet-enriched clusters wholesale (vs. the
house flag-don't-drop), and `emptyDrops` does empty-droplet calling but not ambient-RNA
denoising — use CellBender/decontX if you need decontamination.
