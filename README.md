# sc-skills

An opinionated, **Seurat-first** core pipeline of [Claude Code](https://www.claude.com/product/claude-code) *Agent Skills* for single-cell RNA-seq analysis — from raw CellRanger output through integration and marker-verified annotation to publication-grade figures, under a shared reproducibility contract.

Skills are plain `SKILL.md` folders, so they load in Claude Code and (unmodified) in Codex and Cursor skill directories.

## What's inside (7 skills)

**Orchestration**
- **single-cell** — entry point; orients an analysis and routes to the right stage.

**Core pipeline**
- **sc-preprocessing** — CellRanger output → analysis-ready Seurat: ambient RNA (CellBender / decontX), doublet flagging, MAD-based QC, gene-symbol standardization.
- **sc-integration** — batch correction and integrated embeddings; Harmony vs scVI/scANVI, with integration-quality checks.
- **sc-annotation** — cell-type labels via reference transfer, pretrained models, SingleR, or LLM marker-prompting — every label verified against canonical markers.

**Figures & reproducibility standards**
- **scientific-plotting** — R-first publication plotting (tidyplots, scop, ggplot2 + cowplot/ggh4x, pheatmap/ComplexHeatmap).
- **sc-conventions** — single-cell house rules (Seurat as source of truth, the Seurat↔AnnData bridge, palettes, cell-type ordering, figure sizing), layered on top of…
- **scientific-reproducibility** — the universal output contract: dual-version figures, source-data export, exact p-values (not stars), publication dpi/vector, a running `methods.md` log, and parameter-encoded output filenames.

This is a deliberately minimal backbone — the pipeline stages and standards that every single-cell analysis hangs off. Downstream methods (trajectory, gene regulatory networks, cell–cell communication, spatial, CRISPR screens, differential abundance, pseudobulk DE) are intentionally out of scope here.

## Install

**As a Claude Code plugin (recommended):**

```bash
claude plugin marketplace add mianaz/sc-skills
claude plugin install sc-skills@sc-skills
```

**Or copy the skill folders manually:**

```bash
git clone https://github.com/mianaz/sc-skills.git
cp -R sc-skills/skills/* ~/.claude/skills/
```

Keep all seven together — they cross-reference each other, and `sc-conventions` + `scientific-reproducibility` underlie every stage.

## Notes

Referenced software (Seurat, scanpy, scVI, scDblFinder, harmony, …) are ordinary R/Python packages — install the ones a step uses; they are libraries, not skills. A few figure recipes cite the published papers their style was drawn from (e.g. *Popescu et al. Nature 2019*, *Cao et al. Science 2020*).

## Extending the backbone

This is a starting point, not a ceiling. To add a pipeline stage or a method skill:

1. Drop a `skills/<name>/SKILL.md` folder in (plus any `references/*.md`).
2. Have it **inherit the contract** — point to `scientific-reproducibility` (dual-version figures, source data, exact p-values, a `methods.md` log) and, for single-cell work, `sc-conventions`.
3. Add one line to `single-cell`'s **Routing** section so the orchestrator can find it.
4. Keep the activation `description` bounded to a single intent, so it doesn't collide with an existing skill.

Downstream methods — trajectory, gene regulatory networks, cell–cell communication, spatial, CRISPR screens, differential abundance, pseudobulk DE — all follow this pattern and are the natural next additions.

## License

MIT © 2026 Miana ([@mianaz](https://github.com/mianaz)).

Recipes drawn from published figures cite their sources; those methods remain the intellectual property of their respective authors, referenced here under fair-use scholarly attribution.
