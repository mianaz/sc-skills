# Visium hex-grid cell-type colocalization

For each Visium spot, compute a weighted colocalization score between two cell types (e.g. B and T cells) using the spot itself and its 6 hexagonal neighbors. Threshold candidate spots by mean + k×SD on colocalization AND per-type proportion. Remove isolated spots via hierarchical clustering (max-linkage, h=1). Generalizes to any two deconvolved cell-type proportions.

(source: sc-paper-distill/papers/2026-pan-cancer-TLS.md M1)

```r
# Visium 55µm hex grid: 6 neighbors per spot (col ±2 same row, col ±1 row ±1)
get_hex_neighbors <- function(df, spot_idx) {
  x <- df[spot_idx, "col"]; y <- df[spot_idx, "row"]
  which((df$col == x-2 & df$row == y) | (df$col == x-1 & df$row == y-1) |
        (df$col == x-1 & df$row == y+1) | (df$col == x+2 & df$row == y) |
        (df$col == x+1 & df$row == y-1) | (df$col == x+1 & df$row == y+1))
}

# Compute colocalization: self-product + 0.5× neighbor cross-products
for (m in seq_len(nrow(spot_df))) {
  nbrs <- get_hex_neighbors(spot_df, m)
  score <- spot_df[m, "typeA"] * spot_df[m, "typeB"]
  for (n in nbrs)
    score <- score + 0.5 * spot_df[m, "typeA"] * spot_df[n, "typeB"]
  spot_df$coloc[m] <- score
}

# Threshold: mean + 0.5×SD for colocalization, mean - 0.5×SD for each proportion
candidates <- spot_df[spot_df$coloc >= mean(spot_df$coloc) + 0.5*sd(spot_df$coloc) &
                      spot_df$typeA >= mean(spot_df$typeA) - 0.5*sd(spot_df$typeA) &
                      spot_df$typeB >= mean(spot_df$typeB) - 0.5*sd(spot_df$typeB), ]

# Remove isolated single spots via max-linkage clustering
coor <- candidates[, c("row", "col")]
coor$label <- cutree(hclust(dist(coor, "maximum"), "single"), h = 1)
# Keep only spots with ≥1 hex neighbor in the candidate set
```
