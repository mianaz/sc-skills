# R analysis-script organization

How to structure the `analysis*.R` compute scripts in an analysis category. The
two-script default is `analysis.R` (compute) + `plotting.R` (figures from results/).
These rules say when and how to grow beyond it.

## De-stage with checkpoints, not env switches

Prefer a **linear, end-to-end** script that a scientist can read top to bottom and
run in one go. Make each expensive stage **resumable** with a file-existence
checkpoint rather than gating the whole script behind an env switch (`STEP=`,
`RUN_SECTIONS=`):

```r
ckpt <- file.path(SRC, "denovo_all_markers.csv")
if (file.exists(ckpt)) {
  message("== MARKERS: skip (checkpoint exists) ==")
  mk <- read.csv(ckpt)
} else {
  message("== MARKERS: compute ==")
  mk <- FindAllMarkers(o, ...)
  write.csv(mk, ckpt, row.names = FALSE)
}
```

Re-running then recomputes only what's missing; deleting one checkpoint recomputes
just that stage. A cheap, idempotent final step (e.g. applying a manual label map)
can run **every** time, so editing it + re-running takes effect without recompute.
Reference implementation: `liver_organoid/analysis/atlas_annotation/analysis.R`.

Env switches are still fine for a single binary mode an outer job sets once
(`SCRNA_MODE=full`), but not for slicing one analysis into N mutually exclusive
sections — that's what checkpoints + splitting are for.

## When to split, and how many

Split an analysis script when its **main body exceeds ~500 lines**, counting only
real code — exclude comments, blank lines, `message()`/`cat()` logging, and
`library()`/env/dir setup. Up to **5** coherent scripts per category when the
analysis genuinely uses multiple methods/phases; otherwise keep fewer.

Name split scripts `analysis_<phase>.R` by what they DO, in run order, e.g.
`analysis_export.R` → `analysis_transfer.R` → `analysis_finalize.R` →
`analysis_diagnostics.R`. State the run order in each header and the README.
`plotting.R` and a single `.py` (for steps R can't do) sit alongside and don't
count against the 5.

## Comments live in methods/README, not in the script

Strip decorative banner comments (`# =====`) and multi-paragraph rationale walls.
A script gets a short header (purpose, phase/run-order, inputs, outputs) and brief
inline notes only for genuinely non-obvious steps (a fragile hack, a non-default
parameter and why). The scientific rationale, parameter justifications, and
per-result narrative belong in `methods.md` and the category `README.md`.

## Preserve logic when refactoring

Restructuring (de-staging, splitting, de-commenting) must not change results: keep
the same function calls, parameters, and output paths. Verify by re-running — every
checkpoint should skip and outputs should be byte-stable.
