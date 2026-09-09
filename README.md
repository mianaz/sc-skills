# sc-skills

An opinionated, **Seurat-first** set of [Claude Code](https://www.claude.com/product/claude-code) *Agent Skills* for single-cell RNA-seq analysis — from finding a public accession or raw CellRanger output through integration, marker-verified annotation, condition tests, and common downstream methods, under a shared reproducibility contract — plus a **paper-distillation engine** that grows the suite from the literature.

Skills are plain `SKILL.md` folders, so they load in Claude Code and (unmodified) in Codex and Cursor skill directories.

## What's inside (18 skills)

**Orchestration**
- **single-cell** — entry point; orients an analysis and routes to the right stage.
- **sc-conventions** — single-cell house rules (Seurat as source of truth, Seurat↔AnnData bridge, palettes, cell-type ordering, figure sizing), layered on **scientific-reproducibility**.
- **scientific-reproducibility** — universal output contract: dual-version figures, source-data export, exact p-values, publication dpi/vector, a running `methods.md` log, parameter-encoded filenames.

**Ingress & core pipeline**
- **omic-catalog** — find, browse, and download public omics datasets (GEO, SRA, CellxGene).
- **sc-dataretrieval** — validate a known accession, pick a deposited format tier, load into Seurat.
- **sc-preprocessing** — CellRanger (or retrieved) counts → analysis-ready Seurat: ambient RNA (CellBender / decontX), doublet flagging, MAD-based QC. Per-sample; flag, don't drop.
- **sc-integration** — Harmony, scVI/scANVI, or scGPT embeddings, carried back into the Seurat object, with integration-quality checks.
- **sc-annotation** — candidate labels via reference transfer, pretrained models, SingleR, or LLM marker-prompting — every label verified against canonical markers.

**Condition tests** (keep split: expression vs composition)
- **sc-pseudobulk** — sample-level DE and sample PCA (cells are measurements, not replicates).
- **sc-differential-abundance** — cell-type or neighbourhood proportion shifts (Milo, sccomp, propeller).

**Downstream methods**
- **sc-trajectory** — pseudotime, RNA velocity, cell fate.
- **sc-cellchat** — ligand–receptor communication (CellChat).
- **sc-grn** — TF regulons (pySCENIC, SCENIC+, AUCell).
- **sc-spatial** — Visium / MERFISH / Xenium / related spatial assays.
- **sc-crispr** — Perturb-seq / CROP-seq / Mixscape.
- **sc-target** — scRNA + GWAS / scDRS / Open Targets.

**Figures**
- **scientific-plotting** — R-first publication plotting (tidyplots, scop, ggplot2 + cowplot/ggh4x, pheatmap/ComplexHeatmap).

**Growth engine**
- **sc-paper-distill** — turn a paper + its code into a structured, provenance-tracked digest ending in *Promotion Proposals* that fold vetted recipes back into the suite. User-invoked (call it by name). Ships with several worked example digests; the local marker/paper CSV databases stay gitignored.

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

Keep the suite together — skills cross-reference each other, and `sc-conventions` + `scientific-reproducibility` underlie every stage.

## Notes

Referenced software (Seurat, scanpy, scVI, scDblFinder, harmony, CellChat, …) are ordinary R/Python packages — install the ones a step uses; they are libraries, not skills. Figure and method recipes cite the published papers they were drawn from.

## Extending the suite

**Automated — the distillation loop.** Hand `sc-paper-distill` a paper + its code repo. It writes a digest and proposes promotions. You approve the ones worth keeping; they land in `scientific-plotting`, `sc-conventions`, or a method skill with a back-link to the source paper.

**Manual — add a skill.** Drop `skills/<name>/SKILL.md` (plus `references/` as needed), inherit `scientific-reproducibility` / `sc-conventions`, add one routing line to `single-cell`, and keep the activation `description` to a single intent.

## License

MIT © 2026 Miana ([@mianaz](https://github.com/mianaz)).

Recipes drawn from published figures cite their sources; those methods remain the intellectual property of their respective authors, referenced here under fair-use scholarly attribution.
