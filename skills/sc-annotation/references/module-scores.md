# Gene-set scores with UCell

Apply the [shared source and reuse rules](../../sc-conventions/references/gene-sets.md)
before running this recipe. Prefer `UCell::AddModuleScore_UCell` for per-cell
signatures. It ranks expression within each cell, so the unsmoothed score for a
fixed cell, gene universe and signature does not depend on which other cells are
included. See [Andreatta and Carmona, 2021](https://doi.org/10.1016/j.csbj.2021.06.043)
and the [current UCell manual](https://bioconductor.org/packages/release/bioc/manuals/UCell/man/UCell.pdf).

## Score a saved set

The example assumes a complete, unweighted, species-matched record selected from
the shared library. `signature_file` is its path. Use the full RNA gene universe
with raw counts or nonnegative log-normalized expression. For a Seurat v5 object
with split counts layers, join the RNA layers before this example. Avoid scaled
or integrated expression and matrices restricted to signature genes or HVGs.

```r
record <- jsonlite::read_json(signature_file, simplifyVector = TRUE)
stopifnot(record$status == "complete", length(record$weights) == 0L)
up <- unique(unlist(record$genes_up))
down <- unique(unlist(record$genes_down))
members <- data.frame(
  original_gene = c(up, down),
  direction = rep(c("up", "down"), c(length(up), length(down))))
members$mapped_gene <- members$original_gene  # same-species example
members$present <- members$mapped_gene %in% rownames(seu[["RNA"]])
members$used_gene <- ifelse(members$present, members$mapped_gene, NA_character_)
stopifnot(any(members$present))

out <- file.path("source_data", "gene_sets", record$set_id)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
file.copy(signature_file, file.path(out, "definition.json"), overwrite = TRUE)
write.csv(members, file.path(out, "genes_original_and_used.csv"), row.names = FALSE)

features <- setNames(list(c(up, paste0(down, "-", recycle0 = TRUE))), record$set_id)
seu <- UCell::AddModuleScore_UCell(
  seu, features = features, assay = "RNA", slot = "counts",
  maxRank = 1500, missing_genes = "impute", w_neg = 1,
  ties.method = "average", BPPARAM = BiocParallel::SerialParam(), name = "_UCell")
score_column <- paste0(record$set_id, "_UCell")
write.csv(data.frame(cell = colnames(seu), score = seu[[score_column]][, 1]),
          file.path(out, "scores.csv"), row.names = FALSE)
saveRDS(list(method = "UCell::AddModuleScore_UCell",
             UCell_version = as.character(packageVersion("UCell")),
             assay = "RNA", layer = "counts", normalization = "raw counts",
             maxRank = 1500, missing_genes = "impute", w_neg = 1,
             ties.method = "average", parallel_backend = "SerialParam", smoothing = "none",
             n_original = nrow(members), n_present = sum(members$present)),
        file.path(out, "scoring_parameters.rds"))
```

`missing_genes = "impute"` retains the source definition and treats absent genes
as zero expression. Inspect coverage before interpreting scores. Choose `maxRank`
using the number of detected features and assay depth; 1500 is the whole-transcriptome
default. Targeted spatial panels need a panel-appropriate value and adequate set
coverage. Save any change in the scoring record.

UCell supports signed sets through `+`/`-` suffixes. Preserve the source directions
and record `w_neg`. Weighted signatures retain their published model unless an
unweighted adaptation is explicitly justified and labeled. Keep unsmoothed scores
as the primary output; label and record any optional neighbor smoothing separately.

## Compare and reuse

Summarize scores within each sample and cell type. Use the biological replicate
for condition tests and record any display scaling separately. Select a plot for
the comparison through `scientific-plotting`.

For a bulk-DE signature projected onto single cells, retain the full bulk DE
table and selection rule in its library record. When matched contrasts exist,
compare its direction with sample-level pseudobulk DE through `sc-pseudobulk`.
Label the single-cell UCell result as an adaptation when the original paper used
another scoring method; retain the original method in the definition record.
