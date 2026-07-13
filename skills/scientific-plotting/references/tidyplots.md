# tidyplots

Use tidyplots as the default for general scientific plots when the data is a tidy dataframe. It is especially good for compact publication-style plots with observations, summaries, exact p-values, split panels, and consistent color adjustments.

Installed local version observed: `tidyplots` 0.4.0. Function families are visible in `help(package = "tidyplots")` and the official reference index.

## Function map

| Need | Functions |
|---|---|
| Raw points | `add_data_points()`, `add_data_points_jitter()`, `add_data_points_beeswarm()` |
| Box/violin | `add_boxplot()`, `add_violin()` |
| Mean/median/count/sum summaries | `add_mean_*()`, `add_median_*()`, `add_count_*()`, `add_sum_*()` |
| Error bars/ribbons | `add_sem_errorbar()`, `add_ci95_errorbar()`, `add_sd_errorbar()`, `add_sem_ribbon()`, `add_ci95_ribbon()` |
| Time course | `add_mean_line(group = ...)`, `add_mean_dot()`, `add_sem_ribbon()` |
| Composition | `add_barstack_absolute()`, `add_barstack_relative()`, `add_areastack_absolute()`, `add_areastack_relative()` |
| Heatmap | `add_heatmap(scale = "none" / "row" / "column")` |
| Labels | `add_data_labels()`, `add_data_labels_repel()` |
| Statistics | `add_test_pvalue()`, `add_test_asterisks()` |
| Axes/size/facets | `adjust_x_axis()`, `adjust_y_axis()`, `adjust_size()`, `split_plot()` |
| Colors | `adjust_colors()`, `colors_discrete_*`, `colors_continuous_*`, `colors_diverging_*` |
| Save | `save_plot()` |

## Observations plus mean/SEM plus exact p-value

```r
library(dplyr)
library(tidyplots)

p <- df |>
  mutate(group = factor(group, levels = c("Control", "Treatment"))) |>
  tidyplot(x = group, y = value, color = group, width = 45, height = 40) |>
  add_data_points_beeswarm(shape = 1, cex = 2.5, corral = "wrap") |>
  add_mean_bar(alpha = 0.35) |>
  add_sem_errorbar() |>
  add_test_pvalue(
    method = "wilcox_test",
    p.adjust.method = "BH",
    label = "p = {format_p_value(p.adj, 0.0001)}",
    padding_top = 0.18,
    step.increase = 0.12,
    hide_info = TRUE
  ) |>
  adjust_x_axis(title = NULL) |>
  adjust_y_axis(title = "Signal (a.u.)", padding = c(0.04, 0.16)) |>
  adjust_legend_position("none")
```

Use `method = "t_test"` only when the design and assumptions justify it. Use `paired_by = subject_id` for paired tests.

## Box or violin plot with points

```r
p <- df |>
  tidyplot(x = condition, y = score, color = condition, width = 55, height = 42) |>
  add_violin(alpha = 0.25, trim = FALSE) |>
  add_boxplot(alpha = 0.25, show_outliers = FALSE, box_width = 0.25) |>
  add_data_points_beeswarm(size = 0.8, cex = 2.5, rasterize = nrow(df) > 1000) |>
  adjust_x_axis(labels = function(x) stringr::str_wrap(x, 16), rotate_labels = 45) |>
  adjust_y_axis(title = "Score") |>
  adjust_legend_position("none")
```

## Paired points

```r
p <- df |>
  tidyplot(x = visit, y = value, color = group) |>
  reorder_x_axis_levels("Baseline", "Week 4", "Week 12") |>
  add_line(group = subject_id, color = "grey75", linewidth = 0.25) |>
  add_data_points(size = 1) |>
  add_test_pvalue(method = "wilcox_test", paired_by = subject_id, hide_info = TRUE)
```

## Time course with mean and SEM ribbon

```r
p <- df |>
  tidyplot(x = day, y = response, color = treatment, width = 70, height = 42) |>
  add_mean_line(group = treatment, linewidth = 0.5) |>
  add_mean_dot(size = 1.5) |>
  add_sem_ribbon(alpha = 0.18) |>
  adjust_x_axis(title = "Day") |>
  adjust_y_axis(title = "Response")
```

## Composition, stacked bars, and stacked areas

```r
# Relative composition by sample or condition
p_bar <- comp |>
  tidyplot(x = sample, y = fraction, color = cell_type, width = 90, height = 45) |>
  add_barstack_relative() |>
  adjust_x_axis(labels = function(x) stringr::str_wrap(x, 12), rotate_labels = 45) |>
  adjust_y_axis(title = "Fraction", labels = scales::percent)

# Time-varying composition
p_area <- comp_time |>
  tidyplot(x = timepoint, y = fraction, color = cell_type) |>
  add_areastack_relative()
```

If proportions were already computed and must not be re-summarized, verify the plotted y variable sums to 1 within each x group.

## Small heatmap

Use tidyplots heatmaps for compact tidy tables, not for complex annotated matrices.

```r
p <- heat_df |>
  tidyplot(x = condition, y = feature, color = z_score, width = 55, height = 70) |>
  add_heatmap(scale = "none", rotate_labels = 45, rasterize = TRUE) |>
  adjust_colors(colors_diverging_BuRd) |>
  adjust_x_axis(labels = function(x) stringr::str_wrap(x, 16))
```

For row/column annotations or split blocks, use `references/complexheatmap.md`.

## Split panels

```r
p <- df |>
  tidyplot(x = condition, y = expression, color = condition, width = 32, height = 28) |>
  add_data_points_beeswarm(size = 0.6, rasterize = nrow(df) > 1500) |>
  add_mean_dash() |>
  add_sem_errorbar() |>
  split_plot(by = gene, ncol = 4, nrow = 2, axes = "margins", axis.titles = "single")
```

Use split panels instead of making one unreadable plot with many genes/categories.

## Saving

```r
readr::write_csv(source_df, "source_data/fig2b_points.csv")

p |> save_plot("figures/fig2b_clean.pdf", width = 90, height = 70, units = "mm")
ggplot2::ggsave("figures/fig2b_clean.png", p, width = 90, height = 70, units = "mm", dpi = 600)
```

Use `save_plot()` for tidyplots-native PDF/multipage workflows and `ggsave()` when you need explicit PNG/TIFF device settings.
