# Pseudobulk: aggregation → sample PCA → per-cluster DE (R-first)

Seurat object is the source of truth. Aggregate raw counts to the biological replicate, then use bulk machinery. Follows sc-conventions for outputs, source data, and methods.md logging.

## Table of contents
- [1. Aggregate to pseudobulk](#1-aggregate-to-pseudobulk)
- [2. Sample-level PCA / MDS](#2-sample-level-pca--mds)
- [3. Per-cluster DE across conditions (edgeR-QLF)](#3-per-cluster-de-across-conditions-edger-qlf)
- [4. DESeq2 alternative](#4-deseq2-alternative)
- [5. Thresholds & pitfalls](#5-thresholds--pitfalls)

## 1. Aggregate to pseudobulk

Sum **raw counts** per sample (global) or per sample × cluster (cell-type-resolved). Never average normalised values — they don't re-sum to a valid library.

```r
library(Seurat); library(Matrix)
# seu: annotated, integrated object with meta.data$sample, $condition, $celltype
DefaultAssay(seu) <- "RNA"

# Seurat >=5 one-liner (returns a summed-counts assay object)
pb <- AggregateExpression(
  seu, assays = "RNA", slot = "counts",
  group.by = c("sample", "celltype"), return.seurat = FALSE
)$RNA                                   # genes x (sample_celltype) matrix

# Manual equivalent (portable, explicit) — genes x group sparse sum
grp <- interaction(seu$sample, seu$celltype, drop = TRUE, sep = "__")
mm  <- sparse.model.matrix(~ 0 + grp)
colnames(mm) <- sub("^grp", "", colnames(mm))
pb  <- GetAssayData(seu, slot = "counts") %*% mm     # genes x group

# per-group cell counts — used to drop thin pseudobulks
ncells <- as.integer(table(grp))[match(colnames(pb), levels(grp))]
names(ncells) <- colnames(pb)
```

Build a column-level metadata table (sample, condition, celltype, n_cells, plus batch/sex/donor covariates) aligned to `colnames(pb)` — this is the design table for both PCA and DE.

## 2. Sample-level PCA / MDS

For "how do samples group?" (atlas_pseudobulk_pca). Use logCPM on the aggregated matrix, top variable genes.

```r
library(edgeR)
keep <- ncells >= 10                               # drop thin pseudobulks first
y  <- DGEList(pb[, keep]); y <- calcNormFactors(y)
lcpm <- cpm(y, log = TRUE, prior.count = 2)
v  <- head(order(matrixStats::rowVars(lcpm), decreasing = TRUE), 2000)
pc <- prcomp(t(lcpm[v, ]), scale. = TRUE)
# ggplot PC1/PC2 colored by condition, shaped by batch — style via scientific-plotting.
# Also: plotMDS(y) for the edgeR-native distance view.
```

Save the pseudobulk matrix + the PC coordinates as source data (`.csv` + `.rds`) per sc-conventions.

## 3. Per-cluster DE across conditions (edgeR-QLF)

Loop clusters; within each, samples are the replicates. Keep covariates in the design.

```r
library(edgeR)
run_cluster_de <- function(ct, meta, pb, ncells, min_cells = 10, min_reps = 3) {
  cols <- meta$celltype == ct & ncells[rownames(meta)] >= min_cells
  m <- droplevels(meta[cols, ]); mat <- pb[, rownames(m)]
  if (min(table(m$condition)) < min_reps) return(NULL)   # not enough replicates
  design <- model.matrix(~ batch + condition, data = m)   # covariates first
  y <- DGEList(mat); y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y); y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "conditiondisease")       # last term = contrast
  topTags(qlf, n = Inf)$table                              # logFC, PValue, FDR
}
de <- lapply(levels(factor(meta$celltype)),
             \(ct) run_cluster_de(ct, meta, pb, ncells))
```

Volcano per cluster and a combined logFC heatmap → scientific-plotting (`references/volcano.md`). Log the design formula, contrast, filter, and per-cluster replicate counts in methods.md.

## 4. DESeq2 alternative

Same aggregation; DESeq2 wants integer counts and the aligned colData.

```r
library(DESeq2)
sub <- meta$celltype == ct & ncells[rownames(meta)] >= 10
dds <- DESeqDataSetFromMatrix(round(pb[, rownames(meta)[sub]]),
                              colData = meta[sub, ], design = ~ batch + condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "disease", "control"))
```

For the Python route (scanpy → DESeq2), use `pydeseq2` on the aggregated matrix.

## 5. Thresholds & pitfalls

- **≥3 replicates per group per cluster.** Below that, don't report DE p-values — say so.
- **≥10 cells per sample × cluster** before a pseudobulk column is trusted (raise for rare types).
- **Sum counts, not means of logcounts.** The single most common error.
- **Don't reuse cell-level DE (`FindMarkers` across conditions) as a replicate-aware test** — that is the pseudoreplication trap pseudobulk exists to fix (Squair et al. 2021).
- **Confounding:** if every disease sample is also one batch, no method can separate them — check the design table before modelling.
- Pseudobulk answers *expression* shifts; for *proportion* shifts across conditions use a composition/differential-abundance test (Milo/sccomp/propeller), not this.
