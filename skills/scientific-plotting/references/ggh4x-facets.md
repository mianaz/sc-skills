# ggh4x — faceted dotplots & GSEA dotplots with colored strips

ggh4x adds per-facet strip styling (`strip_themed`) and free/spaced facets (`facet_wrap2`,
`facet_grid2`) that base ggplot can't do — used for marker dotplots grouped into colored
cell-type boxes and for GSEA dotplots faceted by gene-set source.

## Marker dotplot faceted by cell type (colored strips)
Build the dotplot data (from `Seurat::DotPlot(...)$data`, or aggregate yourself), map each gene to a
cell-type box, then facet with colored strips:
```r
library(ggplot2); library(ggh4x)
cols  <- ggpubr::get_palette("Paired", k = nlevels(dp$cluster))
strip <- strip_themed(background_x = elem_list_rect(fill = cols),
                      text_x       = elem_list_text(color = "white"))
ggplot(dp, aes(features.plot, id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  facet_wrap2(~ cluster, scales = "free_x", strip = strip, nrow = 1) +
  scale_color_gradient(low = "#ffffff", high = "firebrick3", name = "avg.exp") +
  scale_size(name = "% exp") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank())
```

## GSEA dotplot — pathways × conditions, size=-log10(padj), color=NES
```r
strip_cols <- c(GOBP = "#E74C3C", GOCC = "#3498DB", GOMF = "#2ECC71")
strip <- strip_themed(background_y = elem_list_rect(fill = strip_cols[levels(df$source)]),
                      text_y       = elem_list_text(color = "white", face = "bold", angle = 90))
nes_lim <- max(2, ceiling(max(abs(df$NES), na.rm = TRUE)))
ggplot(df, aes(condition, pathway_clean)) +
  geom_point(aes(size = -log10(padj), color = NES)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        limits = c(-nes_lim, nes_lim), name = "NES") +
  scale_size_continuous(range = c(2, 8), name = expression(-log[10]~P[adj])) +
  facet_grid2(source ~ ., scales = "free_y", space = "free_y", strip = strip) +
  theme_cowplot() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
```

## Pretty pathway names (MSigDB → Title Case)
```r
pretty_pathway <- function(x) x |>
  stringr::str_remove("^(GOBP|GOCC|GOMF)_") |>
  stringr::str_replace_all("_+", " ") |> stringr::str_to_title() |>
  stringr::str_replace_all(c("\\bMhc\\b"="MHC","\\bNk\\b"="NK","\\bDna\\b"="DNA",
                             "\\bCd(\\d+)\\b"="CD\\1"))
# then stringr::str_wrap(..., width = 40) for axis fit
```
Save dual-version + source data (sc-conventions).
