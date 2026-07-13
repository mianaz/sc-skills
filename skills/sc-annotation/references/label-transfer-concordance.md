# Cross-reference label-transfer concordance (Cohen's κ matrix)

When several reference atlases (or annotation methods) can label the same cells, quantify how much they **agree** instead of picking one blindly. Build an N×N Cohen's κ matrix and read its structure: high mutual κ = robust consensus; a reference that disagrees with all others is a candidate outlier (wrong granularity, batch/domain gap, or mislabeled reference).

## 1. Harmonize label vocabularies first
κ is meaningless if references use different names for the same type. Map every reference's labels to **one common (usually coarse) vocabulary** before comparing.

```r
harmonize <- function(x, map) factor(unname(map[as.character(x)]), levels = sort(unique(map)))
# map: named vector, e.g. c("Hepatocyte"="Hep","Hepatocytes"="Hep","KC"="Macrophage",...)
labs <- lapply(pred_list, harmonize, map = common_map)   # pred_list: one label vector per reference
```

## 2. Pairwise κ on shared cells
Transfer is **directional** (reference i → query), so entry (i, j) is not forced to equal (j, i) when each was transferred onto a different target — keep the matrix asymmetric and report it as such (the asymmetry is informative: which reference generalises onto which).

```r
library(irr)   # or psych::cohen.kappa
refs <- names(labs); K <- matrix(NA, length(refs), length(refs), dimnames = list(refs, refs))
for (i in refs) for (j in refs) {
  ok <- !is.na(labs[[i]]) & !is.na(labs[[j]])          # cells labelled by both
  if (sum(ok) > 20)
    K[i, j] <- irr::kappa2(cbind(as.character(labs[[i]][ok]),
                                 as.character(labs[[j]][ok])))$value
}
```

## 3. Read the matrix
- Symmetrize for a consensus view (`(K + t(K))/2`) OR keep directional to see transfer asymmetry.
- Cluster/annotate a heatmap (pheatmap / ComplexHeatmap via scientific-plotting); RdBu or viridis, values in cells.
- **Interpretation:** references that agree strongly with each other but poorly with one outlier suggest the outlier is at a different granularity or from a mismatched domain (e.g. an adult reference transferred onto fetal cells). Low κ everywhere = the common vocabulary is too fine or the cells are genuinely ambiguous → coarsen and re-check.
- κ guide: <0.2 poor, 0.2–0.4 fair, 0.4–0.6 moderate, 0.6–0.8 substantial, >0.8 near-identical.

## 4. Report
Save the κ matrix as source data (`.csv`) and the heatmap dual-version (sc-conventions). Record the common vocabulary map, the shared-cell counts per pair, and which references were concordant vs outliers in methods.md. Consensus labels should lean on the mutually-concordant references; treat an outlier reference's unique calls as hypotheses, not ground truth.
