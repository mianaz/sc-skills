# Single-cell project layout

Inherits the standard `figures/ source_data/ methods.md` layout from **scientific-reproducibility** (dual-version figures + per-plot csv/rds source data live there). Single-cell adds the `objects/` directory:

```
project/
  figures/        # dual-version figures — see scientific-reproducibility
  source_data/    # per-plot .csv + .rds — see scientific-reproducibility
  objects/
    *.rds         # Seurat objects (the source of truth)
    *.h5ad        # intermediate AnnData for Python detours
  methods.md      # running analysis log — see scientific-reproducibility
```

Use a single `<plot>` basename for a figure and its source data so they are trivially linked.
