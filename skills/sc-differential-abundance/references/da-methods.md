# Differential abundance — Milo, sccomp, propeller (R-first)

Seurat/SCE object is the source of truth; n = samples per group. Follows sc-conventions for outputs and methods.md.

## Table of contents
- [1. Milo — neighbourhood DA (cluster-free)](#1-milo--neighbourhood-da-cluster-free)
- [2. sccomp — Bayesian compositional (cluster-level)](#2-sccomp--bayesian-compositional-cluster-level)
- [3. propeller — fast frequentist (cluster-level)](#3-propeller--fast-frequentist-cluster-level)
- [4. Reporting](#4-reporting)

## 1. Milo — neighbourhood DA (cluster-free)

Tests abundance in overlapping KNN neighbourhoods, so it localises *where* on the manifold composition shifts without committing to clusters. NB-GLM (edgeR) per neighbourhood + spatial FDR.

```r
library(miloR); library(SingleCellExperiment)
# sce with a batch-corrected reduced dim (e.g. "HARMONY" or "SCVI") + colData$sample, $condition
milo <- Milo(sce)
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "HARMONY")
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE)
plotNhoodSizeHist(milo)                       # aim for peak >= 5x n_samples
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample = "sample")

design <- data.frame(colData(milo))[, c("sample", "condition", "batch")]
design <- distinct(design); rownames(design) <- design$sample
milo <- calcNhoodDistance(milo, d = 30, reduced.dim = "HARMONY")
res  <- testNhoods(milo, design = ~ batch + condition, design.df = design)
# res: logFC, SpatialFDR per neighbourhood. Annotate neighbourhoods to cell types:
res  <- annotateNhoods(milo, res, coldata_col = "celltype")
```
Visualise with `plotNhoodGraphDA()` (DA graph) and a beeswarm of logFC by cell type (`plotDAbeeswarm`) → style via scientific-plotting. `k`/`prop` control neighbourhood size; neighbourhoods should contain many samples or the GLM is unstable.

## 2. sccomp — Bayesian compositional (cluster-level)

Models per-sample cell-type counts as compositional with variability + outlier handling; returns effect sizes with credible intervals and FDR.

```r
library(sccomp)
res <- single_cell_data |>                    # Seurat or SCE or a cell-level data frame
  sccomp_estimate(
    formula_composition = ~ condition,        # add covariates: ~ batch + condition
    formula_variability = ~ 1,
    .sample = sample, .cell_group = celltype,
    cores = 4
  ) |>
  sccomp_test()                               # FDR + credible intervals
sccomp::plot_summary(res)                      # boxplot of proportions + significant contrasts
```
Outputs a table of `c_effect` (composition logit effect), credible interval, and `c_FDR` per cell type. Good default: robust to a few aberrant samples that break frequentist tests.

## 3. propeller — fast frequentist (cluster-level)

Transforms proportions (logit or arcsin-sqrt) then limma moderated t/F. Fast, sensible with few replicates.

```r
library(speckle); library(limma)
# transform = "logit" (default) or "asin"; handles 2-group or >2-group (ANOVA-style)
out <- propeller(clusters = seu$celltype, sample = seu$sample, group = seu$condition,
                 transform = "logit")
# out: per-cell-type proportions per group, Tstatistic, P.Value, FDR
```
For a covariate-adjusted design, build the transformed-proportion matrix with `getTransformedProps()` and fit your own `limma` design (`~ batch + condition`).

## 4. Reporting

- State n samples per group and the method + design formula in methods.md.
- Report effect (logFC / composition effect) **and** FDR; never a bare "proportion went up".
- Save the per-type/neighbourhood result table as source data (`.csv`), and the DA figure dual-version, per sc-conventions.
- Interpret jointly: because proportions are compositional, flag that a significant increase in one type shows as decreases elsewhere.
