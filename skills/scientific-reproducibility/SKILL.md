---
name: scientific-reproducibility
description: "Use to apply or audit the universal output & reproducibility contract for data analysis in R or Python — dual-version figures, source-data export, exact p-values (not stars), publication dpi/vector, a running methods.md log, param-encoded filenames, and routing compute-heavy steps to a job scheduler. Single-cell work inherits this through sc-conventions; use this directly for non-single-cell analysis. NOT for application/library/infrastructure code or throwaway exploration."
---

# Scientific reproducibility — the output contract

## Overview
Analysis outputs must be **rebuildable and auditable by someone who wasn't there** — including future-you. Every figure ships in two forms with its source data; every choice lands in a log; every expensive run states itself on disk. These are correctness/reproducibility rules, not style preferences.

## When to use
Any analysis step that produces a figure, saves a results object, or makes an analytical choice (test, threshold, parameter) — in R or Python. **When NOT:** application/library/infrastructure code (that is general code-simplification territory), or genuinely throwaway exploration you will not report.

## Quick reference
| Concern | Rule | Detail |
|---|---|---|
| Layout | `figures/  source_data/  objects/  methods.md` at project root | below |
| Every figure | save `_titled` (title + full stat details) AND `_clean` (no title); PNG ≥600 dpi + vector PDF | references/figure-output-contract.md |
| Every figure | export the tidy data behind it → `source_data/<plot>.csv` + `.rds` | references/figure-output-contract.md |
| Comparisons | run the appropriate test, show the **exact p-value**, not significance stars (unless stars requested) | references/figure-output-contract.md |
| Logging | append method+version, refs, params+rationale, results to `methods.md` after each step | references/methods-log.md |
| Swept runs | encode key params in output filenames; write an auditable run/QC log table | references/methods-log.md |
| Naming | cluster/sample/category names = letters, digits, `_` only | inline |
| Compute-heavy | never run heavy/GPU steps on a login node — submit via the scheduler | inline |

## Compute-heavy steps go through the scheduler
Alignment, model training/integration, large simulations, GPU jobs, and atlas-scale runs belong in a batch job, not an interactive/login shell — login nodes are shared and may lack the CPU features (e.g. AVX-512) the compute nodes have, so code can silently run slow or crash with opaque illegal-instruction errors. Right-size CPU/memory/time and submit through the job scheduler. On a Slurm cluster, use the `hpc-etiquette` skill to write the script.

## Naming hygiene
Fix names at creation time (annotation, sample sheet), not after outputs exist: letters, digits, `_` only. Slashes, spaces, `@`, and other special characters break file paths, classifiers, and regex parsing with vague downstream errors.

## Precedence over general code skills
These reproducibility rules **override** general code-simplification skills for analysis work: do not "simplify away" dual-version plots, source-data exports, exact p-values, dpi/vector output, or the `methods.md` log. Those skills govern glue code and plumbing; they do not trump analysis correctness or reproducibility.

## Relationship to domain suites
Domain conventions **inherit** this contract and add specifics — e.g. `sc-conventions` (single-cell: Seurat source-of-truth, R↔Python bridge, palettes, cell-type ordering). When a domain conventions skill is active, follow it *and* this contract.
