# Gene-symbol standardization (HGNChelper)

Per-sample, right after loading counts and before any merge. Standardizes aliases/outdated
symbols to current official symbols so genes line up across samples and reference gene sets.

```r
library(HGNChelper)

# species = "human" or "mouse"
check <- suppressMessages(checkGeneSymbols(rownames(dat), species = "mouse"))
# multi-mapping suggestions arrive as "A /// B" — keep the first
check$Suggested.Symbol <- sapply(strsplit(check$Suggested.Symbol, " /// "), `[`, 1)

dat   <- dat[!is.na(check$Suggested.Symbol), ]
check <- check[!is.na(check$Suggested.Symbol), ]
rownames(dat) <- check$Suggested.Symbol[rownames(dat) %in% check$x]
```

## Collapse duplicate symbols after correction
Correction can map two old symbols onto one current symbol — sum their counts:
```r
dups <- unique(rownames(dat)[duplicated(rownames(dat))])
for (g in dups) {
  idx <- which(rownames(dat) == g)
  dat[idx[1], ] <- colSums(dat[idx, ])
  dat <- dat[-idx[-1], , drop = FALSE]
}
```

## Cross-species ortholog mapping (e.g. human cell-cycle genes → mouse)
```r
mmus_s <- gprofiler2::gorth(cc.genes.updated.2019$s.genes,
            source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_s <- checkGeneSymbols(mmus_s, species = "mouse")$Suggested.Symbol
```
Log the species and that symbols were standardized in methods.md.
