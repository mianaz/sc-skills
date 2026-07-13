# Atlas gene filtering (multi-study merges)

When combining many studies into one atlas, restrict to genes that are (1) detected in
enough studies and (2) protein-coding, to drop study-specific noise and non-coding clutter.
Apply at/just before the merge, after per-sample symbol standardization (gene-symbols.md).

## Keep genes detected in ≥50% of studies
```r
gene_list <- lapply(studies, function(st) {
  cells <- colnames(obj)[obj$study == st]
  cm    <- Seurat::GetAssayData(obj, layer = "counts")[, cells]
  rownames(cm)[Matrix::rowSums(cm > 0) > 0]
})
gene_counts  <- table(unlist(gene_list))
min_datasets <- ceiling(length(gene_list) * 0.5)
common_genes <- names(gene_counts[gene_counts >= min_datasets])
```

## Restrict to protein-coding (biomaRt)
```r
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # mouse: mmusculus_gene_ensembl
pc <- biomaRt::getBM(attributes = "external_gene_name",
                     filters = "biotype", values = "protein_coding",
                     mart = ensembl)$external_gene_name
genes_keep <- intersect(common_genes, pc[pc != ""])
obj <- subset(obj, features = genes_keep)
```

## Drop tiny samples
```r
n <- table(obj$sample_id)
obj <- subset(obj, subset = sample_id %in% names(n[n >= 100]))
```
## Retain TR/IG genes when immune cells matter

The protein-coding filter above can silently drop T-cell-receptor and immunoglobulin genes
(some annotations don't biotype them as `protein_coding`). For lymphocyte-rich tissues
(tumour immune microenvironment, lymphoid organs) **union the TR/IG biotypes back in** so
T/B identity and Ig signatures survive:
```r
trig <- biomaRt::getBM(attributes = "external_gene_name", filters = "biotype",
                       values = c("TR_V_gene","TR_J_gene","TR_C_gene","TR_D_gene",
                                  "IG_V_gene","IG_J_gene","IG_C_gene","IG_D_gene"),
                       mart = ensembl)$external_gene_name
genes_keep <- union(genes_keep, intersect(rownames(obj), trig))
```
(source: Xue et al. Nature 2022, kept protein-coding ∪ TR/IG.)

## Stress/ribo gene removal before HVG — alternative school, NOT the default

Some atlas pipelines (incl. scPLC) **delete** MT, HSP, ribosomal (RP[SL]), and
dissociation-induced gene sets from the matrix *before* HVG selection, so those programs
can't drive variable genes. **This is a deliberate alternative to the house approach.** The
house default keeps those genes and instead handles them as QC fractions + regression
covariates (`percent.mt`, cell-cycle) at scaling — see qc.md — i.e. flag/regress, don't delete.
Only switch to pre-HVG removal if stress/ribosomal programs are demonstrably driving clusters
after integration *and* regression didn't resolve it; record the swap and rationale in methods.md.

Record the detection fraction, biotype filter, TR/IG retention, and min-cells cutoff in methods.md.
