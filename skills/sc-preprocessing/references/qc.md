# QC metrics — flag, don't drop (this stage)

## Metrics (human shown; mouse patterns in comments)
```r
seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")     # mouse: "^(Mt|mt)-"
seu[["percent.rb"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]")  # mouse: "^(Rpl|Rps)"
seu[["percent.hb"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^HB[^(P)]")# mouse: "^(Hbb|Hba)-"
```

## Data-driven thresholds (MAD-based, preferred over hardcoded numbers)
Center on the median, scale by MAD; common choice is median ± 5·MAD (asymmetric ok).
```r
m  <- median(seu$nCount_RNA);  s  <- mad(seu$nCount_RNA)
mf <- median(seu$nFeature_RNA); sf <- mad(seu$nFeature_RNA)
thr_nCount   <- c(0,   m  + 5 * s)
thr_nFeature <- c(100, mf + 5 * sf)   # hard floor on genes/cell
thr_mt       <- c(0,   25)            # mito ceiling stays a fixed biological cap
```
Record the computed median/MAD and final cutoffs (+ rationale) in methods.md.

## Flag, don't drop
```r
seu$qc_pass <- seu$nCount_RNA   >= thr_nCount[1]   & seu$nCount_RNA   <= thr_nCount[2] &
               seu$nFeature_RNA >= thr_nFeature[1] & seu$nFeature_RNA <= thr_nFeature[2] &
               seu$percent.mt   <= thr_mt[2]
```
Also tabulate pre/post-QC summaries per sample (ncells, median nCount/nFeature/percent.mt,
doublet_rate). Cells are flagged, not deleted, so thresholds can be revisited after integration.

## Cell-cycle scoring + regression (at normalize/scale)
Score before scaling, then regress S/G2M (and percent.mt) out so cycle doesn't dominate clustering.
```r
seu <- Seurat::NormalizeData(seu)
seu <- Seurat::CellCycleScoring(seu,
         s.features   = cc.genes.updated.2019$s.genes,     # mouse: map via gprofiler2::gorth (gene-symbols.md)
         g2m.features = cc.genes.updated.2019$g2m.genes)
seu <- Seurat::ScaleData(seu, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
# SCT alternative: SCTransform(seu, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
```

Visualize via **scientific-plotting** (violin of nFeature/nCount/percent.mt/percent.hb/percent.rb;
scatter nCount vs nFeature colored by qc_pass), dual-version + source data per sc-conventions.
