---
paper: "hepatic-macrophage-niches — Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches"
authors: "Guilliams M et al."
doi: "10.1016/j.cell.2021.12.018"
journal: "Cell"
year: 2022
repo: "https://github.com/guilliottslab/scripts_GuilliamsEtAll_Cell2022 @ 116cc3f3289972fac6f13cb061a86fbcde07ce58"
data_accession: "GSE192742; PMCID: PMC8809252; preprint DOI: 10.1101/2021.10.15.464432"
distilled: 2026-07-07
tags: [human, mouse, macaque, pig, chicken, hamster, zebrafish, liver, macrophage, kupffer-cell, LAM, CITE-seq, scRNA-seq, snRNA-seq, Visium, spatial-proteomics, MICS, Molecular-Cartography, NicheNet, cross-species]
---

# hepatic-macrophage-niches — distillation

## 1. Overview
- **What it is & why notable:** A cross-technology, cross-species liver atlas focused on macrophage niches. The study combines transcriptomics, proteomics, and spatial profiling to identify bona fide Kupffer cells (KCs), define non-KC macrophage niches (bile duct, capsule, central vein), and derive conserved KC programs across species.
- **Data:** healthy and steatotic human liver plus murine liver; additional species: macaque, pig, chicken, hamster, zebrafish. Modalities include scRNA-seq, snRNA-seq, CITE-seq, Visium spatial transcriptomics, Visium highly multiplexed protein, MICS imaging, Molecular Cartography (smFISH), confocal IF, flow cytometry, qPCR, and perturbation/in vivo genetics.
- **Code availability:** Analysis scripts are available at `https://github.com/guilliottslab/scripts_GuilliamsEtAll_Cell2022` (pinned commit above), including FastCAR ambient correction, Seurat/Harmony preprocessing, scVI/TotalVI notebooks, post-vi integration export, and Visium preprocessing/clustering scripts. Data are deposited in GEO (`GSE192742`) and atlas browsing is at `www.livercellatlas.org`.
- **Evidence level:** `full-text` (local PDF), `full-text` (local supplementary xlsx bundle), and `repo-only` (for code-specific implementation details).

## 2. Methods catalog

### M1. Multi-omic atlas assembly with modality-aware latent integration
- **Purpose:** build robust liver cell identity maps despite protocol and modality differences.
- **How used:** CellRanger processing, ambient RNA correction (FastCAR cutoff 0.05), QC by MAD metrics + PCA outlier removal, SCTransform + optional Harmony, then `scVI` (non-CITE) or `TotalVI` (CITE) for denoised expression/protein and clustering.
- **On what data:** pooled mouse and human cell/nuclei datasets across isolation protocols.
- **Tools / packages:** CellRanger 3.1.0, Scater, Seurat, Harmony, scvi-tools (scVI/TotalVI), pheatmap, ggplot2.
- **Code:** `2_FastCAR.R:27,58,128`; `3a_SCT.R:100-120,243,305-306,318`; `3b_Harmony.R:106-126,269,319-320,332`; `4b_TotalVI.ipynb` (model workflow, HVG=4000); `4c_TotalVI_CiteSeq+non-CiteSeq.ipynb` (mixed CITE/non-CITE); `4d_ScVI.ipynb` (RNA-only model workflow); `4e_post_TotalVI-ScVI.R`.
- **Interpretation:** using modality-aware variational models was central for joint RNA+ADT interpretation and for consistent marker extraction.
- **Reusable?** yes → `sc-integration` / `sc-preprocessing` (decision rule: scVI for RNA-only, TotalVI for CITE).

### M2. Cross-technology harmonization strategy (scRNA + snRNA + CITE + spatial)
- **Purpose:** compare cell types across technologies that capture different molecular layers.
- **How used:** derive denoised gene/protein values from scVI/TotalVI; map cell-type signatures onto Visium zonation; transfer zonation between modalities using amortized latent model; convert CITE ADT denoised matrix to FCS for flow-style gating.
- **On what data:** matched liver samples and pooled atlases in mouse/human.
- **Tools / packages:** Seurat spatial workflow, custom probabilistic modeling (PyTorch), flowCore `write.FCS`.
- **Code:** `4a_pre_TotalVI-ScVI.R:36-50,68-107,167-171` (RNA/ADT export and harmonization prep), `4e_post_TotalVI-ScVI.R:27-73` (denoised genes/proteins, UMAP metadata export), `5_visiumData.R:22-30,94-96,179-182` (spatial workflow bridge).
- **Interpretation:** strategy explicitly bridges sequencing, protein imaging, and cytometry space.
- **Reusable?** yes → `sc-spatial` (signature-to-spatial transfer + latent zonation transfer pattern).

### M3. Cross-species KC signature construction
- **Purpose:** define conserved KC identity rather than species-specific markers.
- **How used:** identify human and mouse KC markers separately; filter by off-target expression using scaled denoised heatmaps; map orthologs via HGNC BioMart + MGI; retain markers that are HVG-compatible across species; refine to compact conserved set and project to additional species by ortholog lookup.
- **On what data:** human and mouse myeloid atlases, then macaque/pig/chicken/hamster/zebrafish.
- **Tools / packages:** scVI/TotalVI denoised matrices, heatmap filtering, HGNC BioMart, MGI.
- **Code:** method logic is described in STAR Methods; downstream marker processing scripts include `3a_SCT.R:342` and `3b_Harmony.R:354` (`FindAllMarkers`), with species-wise KC marker supplements in `mmc9`.
- **Interpretation:** filtering strategy avoids over-general markers and emphasizes cross-species portability.
- **Reusable?** yes → new sub-skill candidate (`sc-cross-species-signatures`) or integrate into `sc-annotation`.

### M4. Spatial deconvolution and zonation as probabilistic latent variables
- **Purpose:** infer cell-type composition and zonation from Visium RNA/protein.
- **How used:** negative-binomial generative model with foreground/background components for multiplexed proteins; latent `z ~ Uniform(0,1)` zonation with spline basis (10 knots, sigma 0.05), variational inference, likelihood-ratio significance for cell-type presence in spots.
- **On what data:** human and mouse Visium transcriptomics + highly multiplexed protein.
- **Tools / packages:** custom PyTorch model (ELBO + Adam), conceptually similar to cell2location/scVI family.
- **Code:** Visium preprocessing/clustering and QC proxies are implemented in `5_visiumData.R:122,141-144,161-162,179-182,185-186`; custom probabilistic deconvolution/zonation model from STAR Methods is not present as a standalone script in this repository.
- **Interpretation:** gives transferable zonation and condition-interaction testing (low vs high steatosis).
- **Reusable?** yes → `sc-spatial` method pattern (latent zonation transfer + LR-based presence testing).

### M5. Differential NicheNet with region-specific prioritization
- **Purpose:** infer niche-specific ligand-receptor programs that explain KC vs non-KC niche biology.
- **How used:** compare defined sender/receiver groups across niches (KC vs central-vein mac vs capsule mac vs LAM); include periportal region-specific weighting for ligand prioritization; average mouse and human prioritization into conservation score and visualize top LR pairs (circos).
- **On what data:** integrated sc/snRNA for mouse and human niches.
- **Tools / packages:** Differential NicheNet (`saeyslab/nichenetr` workflow).
- **Code:** paper repo provides upstream inputs and marker/sender definitions (`4a_pre_TotalVI-ScVI.R`, `4e_post_TotalVI-ScVI.R`); Differential NicheNet implementation is referenced from `https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet.md`.
- **Interpretation:** supports claim that ALK1-BMP9/10 axis is a conserved KC niche signal.
- **Reusable?** yes → `sc-cellchat`/`sc-target` adjacent; add a dedicated recipe under communication analysis skilling.

### M6. Functional perturbation validation of predicted niche axis
- **Purpose:** test causality for inferred signaling.
- **How used:** BM monocyte ac-LDL stimulation + qPCR, Acvrl1 conditional models, Clec4f-Dtr chimera setup, Fc-trap perturbations (ALK1Fc, TGFBRIIFc) with flow and microscopy readouts.
- **On what data:** murine in vivo and ex vivo validation assays.
- **Tools / packages:** flow cytometry, qPCR, confocal IF.
- **Code:** not code-centric; experimental protocol.
- **Interpretation:** strengthens computational-inference-to-function chain.
- **Reusable?** no (biological protocol), but good benchmark for evidence-strength criteria in distills.

## 3. Figure catalog

### F1. Figure 1 — Multi-omic atlas assembly and spatial grounding
- **Type:** UMAP + spatial clusters + zonation maps + multiplex images.
- **What it communicates:** atlas quality and modality alignment (RNA, ADT, spatial transcript, spatial protein).
- **Data shown:** pooled sc/snRNA + CITE with Visium and multiplex imaging overlays.
- **Visual style:** panel-dense multi-modal narrative; matched marker coloring between RNA/protein channels.
- **Code:** partial code available: atlas preprocessing and embedding generation (`3a_SCT.R:243,305-306,318`; `3b_Harmony.R:269,319-320,332`; `4e_post_TotalVI-ScVI.R:27-73`) plus Visium plotting scaffolds (`5_visiumData.R:22-30,179-182,222-245`).
- **Reusable template:** use house `scop` + `scientific-plotting` multi-panel layout where each modality has matched marker legend.
- **Interpretation:** establishes technical validity before biological claims.
- **Reusable?** yes → `scientific-plotting` (multi-modal consistency recipe).

### F2. Figure 2/3 — Macrophage subtypes and niche localization
- **Type:** myeloid re-cluster UMAP, DEG/DEP panels, confocal/Molecular Cartography niche images.
- **What it communicates:** KCs and LAMs are distinct and spatially segregated across bile duct/capsule/central vein contexts.
- **Data shown:** transcript/protein markers plus spatial localization.
- **Visual style:** combine unbiased cluster maps with targeted IF/smFISH confirmation.
- **Code:** partial: cluster/marker workflows and plotting scaffolds in `4e_post_TotalVI-ScVI.R:87-109` and `5_visiumData.R:251`; publication figure assembly scripts are not explicitly provided.
- **Reusable template:** cluster -> marker selection -> spatial confirmation triplet.
- **Interpretation:** supports claim of distinct macrophage niches in healthy liver.
- **Reusable?** yes → `sc-annotation` + `sc-spatial` claim-support blueprint.

### F3. Figure 4/S7 — Cross-species bona fide KC identification
- **Type:** cross-species UMAPs, conserved-signature maps, marker validation images.
- **What it communicates:** KC core program is evolutionarily conserved; VSIG4/CD5L-centered identity is robust.
- **Data shown:** human/mouse plus additional vertebrate species profiles.
- **Visual style:** same annotation logic propagated across species panels.
- **Code:** partial: cross-species expression processing scripts are provided; species-level KC marker tables are in supplement `mmc9` (sheets: Macaque/Pig/Hamster/Chicken/Zebrafish KC DEGs, Conserved & Unique KC DEGs); final figure-assembly scripts are not explicitly separated.
- **Reusable template:** conserved-signature transfer with ortholog mapping + per-species sanity plots.
- **Interpretation:** supports portability of KC identity beyond mouse/human.
- **Reusable?** yes → proposed cross-species signature sub-skill.

### F4. Figure 5/6 — Steatosis-driven niche remodeling and stromal context
- **Type:** healthy vs steatotic spatial maps, signature overlays, DEG heatmaps.
- **What it communicates:** LAM localization shifts toward steatotic zones; stromal/vascular context changes with disease.
- **Data shown:** Visium zonation, marker IF, stromal re-clustering.
- **Visual style:** paired condition panels (healthy/steatotic) with same scales.
- **Code:** partial: spatial clustering and differential marker extraction are scripted in `5_visiumData.R:179-182,251`; panel-specific plotting scripts are incomplete.
- **Reusable template:** paired-condition spatial scoring and zonation comparison.
- **Interpretation:** links disease state to niche redistribution.
- **Reusable?** yes → `sc-differential-abundance` + `sc-spatial` integration.

### F5. Figure 6/7/S8 — Differential NicheNet and ALK1-BMP9/10 axis
- **Type:** circos LR plots, feature maps, perturbation validation plots.
- **What it communicates:** conserved ligand-receptor logic predicts and explains KC development/maintenance.
- **Data shown:** LR prioritization plus in vivo/ex vivo validation outputs.
- **Visual style:** computational inference panel adjacent to perturbation readout panel.
- **Code:** analysis recipe available via Differential NicheNet vignette; paper repo includes preprocessing layers needed to generate those inputs (`4a_pre_TotalVI-ScVI.R`, `4e_post_TotalVI-ScVI.R`).
- **Reusable template:** communication inference + perturbation validation evidence ladder.
- **Interpretation:** converts correlational niche maps into mechanistic hypothesis testing.
- **Reusable?** yes → communication-analysis recipe with explicit validation tiers.

## 4. Style observations
- Multi-modality is always shown as a claim chain: cluster map -> marker specificity -> spatial localization -> perturbation/functional test.
- Cross-species comparisons keep a fixed conceptual scaffold (same cell programs, orthologized markers, species-specific confirmation images).
- Condition contrasts (healthy vs steatotic) are mirrored in panel design and scales to reduce interpretability ambiguity.

## 5. Directed outputs requested in this distill

### 5.1 Paper info table

| title | authors | journal | year | doi | main findings summary |
|---|---|---|---:|---|---|
| Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches | Guilliams M et al. | Cell | 2022 | 10.1016/j.cell.2021.12.018 | Built a liver spatial proteogenomic atlas across modalities and species; identified bona fide KC identity and distinct non-KC macrophage niches; showed steatosis-associated niche remodeling; proposed and experimentally supported conserved ALK1-BMP9/10 signaling as a KC niche-development axis. |

### 5.2 Methods synthesis (requested structure)

0) **Study design and key experiments**
- Discovery atlas: scRNA-seq, snRNA-seq, CITE-seq, Visium transcriptomics, Visium highly multiplexed protein, MICS, Molecular Cartography in mouse and human.
- Cross-species extension: macaque, pig, chicken, hamster, zebrafish single-cell/nucleus profiling and ortholog-based signature transfer.
- Functional validation: flow/qPCR, ac-LDL stimulation, conditional Acvrl1 perturbations, chimera and Fc-trap experiments.

1) **Data processing strategy**
- CellRanger for count matrices; FastCAR ambient correction; Scater QC with MAD-based filtering + PCA outlier detection; Seurat preprocessing; Harmony when needed.
- scVI for RNA-only datasets and TotalVI for CITE datasets; denoised outputs used for DEGs/DEPs and marker visualization.
- Visium preprocessing with SCTransform and quality spot filtering; custom probabilistic deconvolution/zonation model.

2) **Cross-technology and cross-species comparison strategy**
- Cross-technology: map cell signatures between RNA and spatial modalities; transfer zonation latent variable between datasets/modalities; convert CITE denoised ADT to FCS for flow-space validation.
- Cross-species: derive stringent human and mouse KC marker sets separately, ortholog-map and intersect/refine, then project to additional species.
- Filtering/integration choices: off-target expression filtering for marker specificity, HVG compatibility check across species, region-specific weighting in Differential NicheNet.

3) **Claim-support analysis and plots**
- Distinct macrophage niches claim: myeloid UMAP + DEG/DEP + confocal/Molecular Cartography localization.
- Conserved KC identity claim: cross-species signature projections, KC marker expression panels, Venn/TF conservation summaries.
- Disease remodeling claim: healthy vs steatotic Visium zonation and LAM/KC signature overlays, NAFLD DE analyses.
- Mechanistic axis claim: Differential NicheNet circos/priority plots + ALK1/BMP9/10 perturbation validation.

4) **Actual reusable code availability**
- **Paper repository (now resolved):** `https://github.com/guilliottslab/scripts_GuilliamsEtAll_Cell2022`.
- **Local assets used in second pass:** full text PDF (`1-s2.0-S0092867421014811-main.pdf`) and supplementary tables (`mmc1`-`mmc10` xlsx) from local folder.
- **Core reusable scripts in repo:**
  - `2_FastCAR.R`: ambient RNA profile estimation and correction workflow.
  - `3a_SCT.R` / `3b_Harmony.R`: Seurat QC, SCTransform/log-normalization, Harmony integration, clustering.
  - `4a_pre_TotalVI-ScVI.R`: RNA/ADT matrix export and CITE/non-CITE harmonization prep.
  - `4b_TotalVI.ipynb`, `4c_TotalVI_CiteSeq+non-CiteSeq.ipynb`, `4d_ScVI.ipynb`: model training/inference templates.
  - `4e_post_TotalVI-ScVI.R`: denoised layers and embeddings extraction for downstream plotting.
  - `5_visiumData.R`: Visium QC, filtering, SCTransform, clustering, spatial plotting, DE extraction.
- **Available external method code used by paper:**
  - Differential NicheNet tutorial workflow: `https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet.md`

### 5.3 Second-pass supplement refinement

- Supplement `mmc9` sheet `Conserved & Unique KC DEGs` contains an explicit 29-gene KC program conserved across chicken/hamster/human/mouse/pig/zebrafish:
  `TXN, CMKLR1, CREG1, ABCA1, C1QA, P2RY13, SLC40A1, FTH1, HMOX1, CTSS, C1QC, CTSD, LIMS1, RBM47, CSF1R, ABCC5, TCN2, RGL1, TNFRSF21, PSAP, SLC23A2, MERTK, LGMN, NR1H3, C1QB, ASAH1, LIPA, GBP2, CCND1`.
- Supplement `mmc2/mmc8/mmc9` provides per-species and per-compartment KC/LAM DEG/DEP tables used to support cross-species and niche claims.

## 6. Marker and signature capture (omics-focused)

The per-marker rows for this paper were appended to persistent database:
- `sc-paper-distill/references/omic_marker_db.csv` (new)

Captured categories for this paper include:
- KC identity and conserved KC program markers
- LAM/bile-duct/capsule/central-vein macrophage context markers
- Region-specific sender markers used in Differential NicheNet prioritization
- LR/axis markers tied to functional validation (ALK1/BMP9/10)

## 7. Promotion Proposals

| # | target skill | kind | summary | source |
|---|---|---|---|---|
| P1 | sc-preprocessing | method | Add explicit modality routing rule: use scVI for RNA-only and TotalVI for CITE; carry denoised matrices into marker/heatmap layers. | M1 |
| P2 | sc-spatial | method | Add latent-zonation transfer recipe (train on one cohort/modality and transfer to another) with spline latent variable and LR presence testing concept. | M4 |
| P3 | sc-annotation | method | Add cross-species marker derivation workflow: species-specific marker discovery -> off-target filtering -> ortholog mapping -> reciprocal HVG compatibility filtering. | M3 |
| P4 | sc-cellchat (or communication skill) | recipe | Add Differential NicheNet niche-comparison recipe with region-specific prioritization factor and conservation scoring across species. | M5 |
| P5 | sc-conventions | requirement | Require claim chains to include at least one orthogonal spatial/protein validation panel for niche-localization claims when data permit. | §4/F2/F5 |
| P6 | sc-paper-distill | fix | Extend distill schema to include mandatory marker/signature curation block that writes to persistent marker DB with experiment/species/context/citation fields. | §6 |
