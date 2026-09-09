---
paper: "DESCARTES fetal atlas — A human cell atlas of fetal gene expression"
authors: "Cao J, O'Day DR, Pliner HA, et al. (Trapnell & Shendure)"
doi: "10.1126/science.aba7721"
journal: "Science"
year: 2020
repo: "https://github.com/JunyueC/sci-RNA-seq3_pipeline @ f8f34970c3a5aebe87d8cce2ed22446974398cee (Zenodo snapshot 10.5281/zenodo.4013713, sci-RNA-seq3_pipeline-master.zip, 2020-09-03); trajectory: https://github.com/cole-trapnell-lab/monocle3"
data_accession: "GEO GSE156793 (expression matrices, supplementary files S1–S7, table S5); dbGaP phs002003.v1.p1 (raw reads); companion ATAC = Domcke et al. Science 2020"
distilled: 2026-07-09
tags: [human, fetal, multi-organ, scRNA-seq, sci-RNA-seq3, atlas, cross-species, trajectory, Monocle3, Garnett]
---

# DESCARTES fetal atlas — distillation

## 1. Overview
- **What it is & why notable:** The gene-expression half of DESCARTES — a single-cell atlas of **~4 million cells across 15 human fetal organs** (72–129 days postconception), generated with **three-level combinatorial indexing (sci-RNA-seq3)**. Notable for (a) scale + breadth on one indexing chemistry, (b) a reusable framework for *quantifying cell-type specificity* (SVM cross-validation F1, Garnett, marker specificity scores), (c) systematic **cross-organ analysis of broadly distributed cell types** (blood, endothelial, epithelial), and (d) straightforward **human↔mouse (MOCA) integration** despite >100 My divergence. Companion to Domcke et al. (chromatin accessibility).
- **Data:** *Homo sapiens*, 28 fetuses, 15 organs (cerebrum, cerebellum, adrenal, lung, kidney, heart, eye, liver, pancreas, intestine, muscle, placenta, spleen, stomach, thymus), 112 samples post-filtering, ~4.06M cells post-filtering. Platform: **sci-RNA-seq3** (mostly nuclei; paraformaldehyde-fixed). Shallow depth (~14k raw reads/cell; median 863 UMI nuclei / 524 genes). 77 main cell types, 657 subtypes, 172 initial per-organ types.
- **Code availability:** The repo (`JunyueC/sci-RNA-seq3_pipeline`, = the Zenodo "demultiplexing script and tutorial") ships **only the upstream processing pipeline** — barcode/UMI attach → STAR align → filter → two-pass UMI dedup → gene count (exon+intron) → sparse matrix + cell/gene annotation. **All downstream analysis and every figure panel are code-absent** (Monocle3/Garnett/Harmony/Seurat/DDRTree referenced by method + tool page `cole-trapnell-lab.github.io/monocle-release/monocle3` only). So: pipeline = `repo-only`; clustering/annotation/trajectory/figures = `methods-only` / `captions-only`.

## 2. Methods catalog

### M1. sci-RNA-seq3 upstream processing → gene count matrix
- **Purpose:** Turn raw sci-RNA-seq3 reads (3-level combinatorial barcodes) into a cell × gene sparse matrix analysis packages can ingest.
- **How used:** RT + ligation barcodes extracted from read1 and **error-corrected to nearest barcode at Hamming/edit distance ≤ 1** (reads with unmatched barcodes discarded); UMI + corrected barcodes appended to read2 name. STAR align read2 (`--outSAMstrandField intronMotif`); filter `q > 30`; **two-pass UMI dedup** (pass 1 exact UMI+tagmentation-site match, pass 2 edit-distance-aware); split BAM by cell barcode; **UMI cutoff = 200** to call a cell (cells below discarded). Parallelized via an R fan-out wrapper over per-sample bash steps.
- **On what data:** every organ/sample from FASTQ; sentinel HEK293T/NIH-3T3 spike-in wells used as per-experiment batch controls.
- **Tools / packages:** Python 2.7, STAR 2.5.2b, samtools 1.4, cutadapt 1.8.3, trim_galore 0.4.1; R (Matrix, tidyverse, data.table); scrublet v0.1 for doublets.
- **Code:** `sci3_main.sh` (full driver); `script_folder/UMI_barcode_attach_gzipped_with_dic.py`, `sciRNAseq_count.py`, `rm_dup_barcode_UMI*.py`, `sam_split.py`.
- **Interpretation:** Standard sci-seq preprocessing; the reusable choices are barcode error-correction radius (≤1), the two-pass dedup, and the UMI≥200 gate.
- **Reusable?** yes → reference recipe if we ever ingest raw sci-RNA-seq3; low priority for our Seurat-first workflow (we normally start from matrices).

### M2. Exon+intron gene counting (nucleus-aware)
- **Purpose:** Recover signal from **nuclear/intronic reads** — essential because most of the atlas is *nuclei*, not whole cells.
- **How used:** Counts are tallied separately for exon and intron features per gene, then **summed** (`combine_exon_intron`) so a gene's UMI count = exon + intron. Reads are categorized (perfect/nearest × intersect/combine × exon/gene) and rolled into `all_exon`, `all_intron`, `all_reads`, `unmatched_rate`; final `UMI_count = all_exon + all_intron`.
- **On what data:** all nuclei; the summary fn returns `df_cell`, `df_gene`, sparse `gene_count` → `sci_summary.RData`.
- **Tools / packages:** R (Matrix, tidyverse).
- **Code:** `gene_count_processing_sciRNAseq.R:1-90` — `sciRNAseq_gene_count_summary()` + `combine_exon_intron()`.
- **Interpretation:** For snRNA-seq, dropping introns throws away most of the signal; exon+intron summation is the correct default (parallels CellRanger `--include-introns`, velocyto spliced+unspliced).
- **Reusable?** yes → sc-preprocessing note: when starting from nuclei, count/keep intronic reads (exon+intron), don't exon-only.

### M3. Cell-type specificity quantification (SVM cross-validation F1)
- **Purpose:** A *quantitative*, label-agnostic confidence score for whether an annotated cluster is a real, separable cell type (vs. over-clustering).
- **How used:** Per organ, subsample ~2000 cells/main type; **5-fold cross-validation with a linear-kernel SVM** on the whole transcriptome; the mean CV **F1 = specificity score**. Null = **permuted cell-type labels** (median F1 ~0.17). Real main types median 0.99; subtypes median 0.77 — both ≫ permuted (P < 2.2e−16, Wilcoxon). Same framework applied to subtype adjudication and to keep/merge decisions.
- **On what data:** all 15 organs, main types and 657 subtypes.
- **Tools / packages:** SVM (linear kernel), Monocle3 for clustering.
- **Code:** not provided (methods-only).
- **Interpretation:** Turns "is this cluster real?" into a number with a permutation null — a rigorous alternative to eyeballing marker separation.
- **Reusable?** yes → sc-annotation: add SVM-CV-F1 + permuted-label null as a cluster-validity / over-clustering check.

### M4. Automated preliminary annotation with Garnett (marker-file classifiers)
- **Purpose:** Reproducible, literature-driven first-pass labels independent of manual clustering, and transferable classifiers.
- **How used:** One **Garnett-for-Monocle3** classifier per organ (v0.2.9), marker files compiled from literature (agnostic to their own clustering); trained on full tissue (except cerebrum, 100k-cell subsample). Applied to external data (adult pancreas) to validate transfer. Manual vs Garnett concordance reported per organ; models posted publicly for reuse.
- **On what data:** all 15 organs; validation on inDrop/scRNA adult pancreas.
- **Tools / packages:** Garnett (Monocle3), `org.Hs.eg.db` 3.10.0.
- **Code:** not provided (methods-only); trained models on descartes.brotmanbaty.org.
- **Interpretation:** Marker-file classifiers give a defensible, portable annotation layer that cross-checks manual labels.
- **Reusable?** yes → sc-annotation: Garnett marker-file classifier as an automated cross-check alongside SCimilarity/Azimuth/SingleR.

### M5. Cross-organ reclustering of broadly distributed compartments
- **Purpose:** Study one cell class (blood / endothelial / epithelial) *across all organs at once* to find organ-specific specializations vs shared programs.
- **How used:** Extract annotated cells of the class from all organs; **downsample to ≤5000 cells/type/organ** to mitigate sampling imbalance; PCA on a gene set built from **top per-type markers** (`q<0.05`, ≥2× fold vs 2nd-ranked, ordered by median q across organs) → UMAP (`max_components=2, n_neighbors=50, min_dist=0.1, metric='cosine'`) → **Louvain (`louvain_res=1e-4`)**. Blood: 103,766 cells; endothelial: 89,219; epithelial: 282,262. For blood, co-embed with an external fetal-liver blood scRNA atlas via **Seurat v3 anchors** (30 dims, top 3000 shared HVGs).
- **On what data:** the three broadly distributed compartments.
- **Tools / packages:** Monocle3 (UMAP/Louvain, `differentialGeneTest()`), Seurat v3 (`FindIntegrationAnchors`/`IntegrateData`), Harmony v1.0 (subclustering batch correction).
- **Code:** not provided (methods-only).
- **Interpretation:** The downsample-then-marker-gene-set recipe is the reusable trick — it prevents the biggest organ from dominating a pan-organ embedding.
- **Reusable?** yes → sc-integration/sc-conventions: pan-compartment cross-tissue embedding recipe (per-type downsample + marker-gene-set PCA).

### M6. Human↔mouse developmental-atlas integration (two ways)
- **Purpose:** Bridge human fetal (72–129d) to mouse organogenesis (MOCA, E9.5–E13.5) despite species + stage gaps.
- **How used:** (a) **Cell-type cross-matching** (Cao et al. 2019 method) maps each of 77 human types to a single MOCA trajectory/subtrajectory; (b) **direct co-embedding**: sample ~100k MOCA + ~65k human (≤1000/type), **Seurat v3 anchors** (30 dims, 3000 shared HVGs), annotate mouse cells by **k=3 human nearest neighbors**. NNLS regression coefficient <0.6 flags no-1:1-match (mouse placenta, skin, gonads).
- **On what data:** all compartments (hematopoietic, endothelial, epithelial done separately too).
- **Tools / packages:** Seurat v3, custom cross-matching, NNLS.
- **Code:** not provided (methods-only).
- **Interpretation:** Ortholog-HVG Seurat anchoring is a "straightforward" cross-species bridge; kNN label transfer + NNLS gate gives a principled no-match call.
- **Reusable?** yes → sc-integration: cross-species anchor + kNN transfer + NNLS no-match gate.

### M7. DDRTree pseudotime for macrophage/microglia branches
- **Purpose:** Order the shared-progenitor → 3 macrophage fates (microglia, phagocytic, perivascular) and find branch-DE TFs.
- **How used:** Top 500 HVGs → 3 PCs → **DDRTree** (`param.gamma=120, norm_method="log", residualModelFormulaStr="~ sm.ns(Total_mRNAs, df=3)"`) → 3 branch trajectories → `kmeans(k=10)`, progenitor = lowest mean development time (central); per-cell pseudotime = distance from progenitor; branch DE via `differentialGeneTest()`.
- **On what data:** 4327 reanalyzed mouse embryonic microglia/macrophages (+ human erythropoiesis trajectory analogously).
- **Tools / packages:** Monocle 2/3 (DDRTree).
- **Code:** not provided (methods-only).
- **Interpretation:** Explicit regression-out of total-mRNA (depth) via natural spline is the copyable parameter choice.
- **Reusable?** partial → sc-trajectory: note the depth-regression `residualModelFormulaStr` idiom; DDRTree is legacy (prefer Monocle3 graph/CellRank now).

## 3. Figure catalog

### F1. Fig 1C — per-organ UMAP grid (15 panels), cluster-labeled in situ
- **Type:** faceted small-multiple UMAP grid (one UMAP per organ), text labels placed on clusters.
- **What it communicates:** the cell-type inventory of each organ at a glance; a "table of contents" for the atlas.
- **Data shown:** each organ's cells in its own 2D UMAP, colored + directly labeled by main cell type.
- **Visual style:** dense grid (4×4-ish), **on-plot text labels next to each cluster** (no side legend), consistent per-organ palette, minimal axes. Summary text block ("15 organs / 112 samples / 4M cells / 77 main cell types") fills the empty grid cell.
- **Code:** caption-only — no code provided.
- **Reusable template:** reconstructed in house style (scop) below (§F1 snippet in proposals).
- **Interpretation:** organs are heterogeneous but share recurrent broadly distributed types (endothelial, blood) — motivating the cross-organ analyses.
- **Reusable?** yes → scientific-plotting: labeled per-group UMAP small-multiple grid recipe.

### F2. Fig 2 — specificity/validity panels (confusion matrices + F1 boxplot + correlation heatmap)
- **Type:** recall confusion matrix (heatmap), F1 box plot (permuted vs main vs subtype), cell-type correlation heatmap (normalized β).
- **What it communicates:** annotations are quantitatively valid and separable, not over-clustered.
- **Data shown:** predicted×actual recall (0–1 viridis); F1 distributions across three conditions; human-subtype × mouse-type correlation β (row-normalized).
- **Visual style:** **viridis** sequential (0→1), square cells, ordered so the diagonal reads; box plot with permuted-null baseline; heatmap row-max normalization.
- **Code:** caption-only — no code provided.
- **Reusable template:** ComplexHeatmap/pheatmap confusion matrix (square tiles, white borders) + tidyplots F1 box with permuted null.
- **Interpretation:** high on-diagonal recall + F1 ≫ permuted ⇒ trustworthy labels.
- **Reusable?** yes → scientific-plotting: cross-validation confusion-matrix + specificity-F1 figure recipe.

### F3. Fig 3B — surface/secreted/TF/ncRNA expression heatmap across 77 main types
- **Type:** 4-panel gene×cell-type expression heatmap (surface proteins | secreted | TFs | ncRNAs).
- **What it communicates:** each main type has a distinct, block-diagonal molecular signature spanning druggable/secreted/regulatory/noncoding gene classes; ncRNAs alone separate developmental groups.
- **Data shown:** rows = 77 main types; columns = top DE genes per class; color = Z-scored, library-size-normalized, log expression **capped to [0,3]**.
- **Visual style:** **block-diagonal ordering** (types ordered so signatures form a staircase), Z-score with an explicit cap, 4 side-by-side panels sharing the row order, sequential fill.
- **Code:** caption-only — no code provided.
- **Reusable template:** ComplexHeatmap staircase (rows ordered by argmax, Z-capped) — reconstructed in house style.
- **Interpretation:** cell identity is legible even restricting to therapeutically/regulatorily interesting gene classes.
- **Reusable?** yes → scientific-plotting: block-diagonal marker heatmap with Z-cap (already partly in suite; strengthen the [0,3] cap + argmax ordering).

### F4. Fig 4D — pan-organ blood marker dot plot
- **Type:** dot plot (genes × blood cell types), two selected markers per type.
- **What it communicates:** newly identified pan-organ markers that specifically label each blood subtype (usable for FACS/labeling), improving on canonical markers that are promiscuous in development.
- **Data shown:** dot **size = % of cells in the type detecting the gene**; dot **color = average expression**; genes grouped/ordered to match cell-type order (staircase).
- **Visual style:** classic Seurat-style dotplot, gene rows ordered to sit on the diagonal of their marking type, sequential color, size legend for detection fraction.
- **Code:** caption-only — no code provided.
- **Reusable template:** scop/ggh4x dotplot pair (already a house recipe) — this is a canonical instance.
- **Interpretation:** e.g. *CD8B*/*CD5* for T, *SORCS1*/*JMY* for ILC3 — cleaner than *CD4*/*CD8A* which leak into macrophages/NK.
- **Reusable?** yes (confirms existing recipe) → scientific-plotting: keep the diagonal-ordered marker dotplot; add these fetal-blood markers as an example set.

### F5. Fig 5A–B — erythropoiesis trajectory + stage marker UMAPs
- **Type:** trajectory-embedded UMAP (Fig 5A) + small-multiple stage-marker UMAPs (Fig 5B) + proportion box plots (5C,I).
- **What it communicates:** a continuous HSPC→erythroid/mega/basophil trajectory; erythroid partitioned into EEP→CEP→ETD stages; unexpected adrenal erythropoiesis.
- **Data shown:** UMAP colored by stage; per-marker UMAPs (expression, Z-scored, aggregated); box+point plots of per-organ blood-cell proportions (samples ≤200 blood cells excluded — an explicit filter).
- **Visual style:** black **directional arrows** overlaid on trajectory UMAP; per-stage marker facet grid; box plots with jittered points and a stated exclusion threshold.
- **Code:** caption-only — no code provided.
- **Reusable template:** scop trajectory UMAP + faceted marker UMAPs; tidyplots box+point with documented exclusion.
- **Interpretation:** conserved erythroid staging (matches mouse Zhu et al.) + adrenal as a minor erythropoietic site.
- **Reusable?** yes → scientific-plotting: trajectory UMAP with directional arrows + faceted stage-marker UMAP recipe.

### F6. Fig 6A–E — human↔mouse co-embedding + conserved-TF marker UMAPs
- **Type:** joint UMAP colored by species / by mouse trajectory / by gestational stage / by human blood type; paired human-vs-mouse expression UMAPs of conserved genes.
- **What it communicates:** human fetal cells project onto matched mouse embryonic trajectories; conserved cell-type-specific TFs mark the same populations in both species.
- **Data shown:** same UMAP recolored 4 ways (species, mouse trajectory, stage, human type); F6E: 2 rows (human/mouse) × 5 types, each colored by a conserved-TF module, **boxed** on the relevant population.
- **Visual style:** **one embedding recolored multiple ways** (small-multiple of the *same* coordinates); grey "other-species" background with colored focal species; **inset boxes** highlighting the matched population; TF names annotated per panel.
- **Code:** caption-only — no code provided.
- **Reusable template:** re-color-the-same-UMAP panel set + grey-background focal overlay (scop `CellDimPlot` split; ggplot grey base + colored subset).
- **Interpretation:** developmental cell-type programs are evolutionarily constrained; mouse is a valid model for these human fetal types.
- **Reusable?** yes → scientific-plotting: "same embedding, recolored N ways" + grey-background focal-overlay recipe.

## 4. Style observations
- **Sequential/viridis** for specificity and expression heatmaps; **Z-score with explicit cap** ([0,3]) is a recurring choice for expression heatmaps and marker UMAPs.
- **Block-diagonal / staircase ordering** of cell types and marker genes so signatures read along the diagonal (heatmaps F3, dotplot F4).
- **On-plot direct labels** instead of side legends for the big per-organ UMAP grid (F1).
- **One embedding, recolored many ways** as the core cross-cutting-analysis idiom (F6, and blood colored by organ vs by type in F4A/B).
- **Grey background = context, color = focal set** (other-species or other-organ cells greyed).
- **Directional arrows** hand-annotate inferred trajectory flow on UMAPs (F4C, F5A).
- **Stated exclusion thresholds in captions** (samples ≤200 blood cells excluded) — provenance discipline worth copying.
- Depth normalization is pervasive: UMI counts library-size-scaled + log + Z before display; total-mRNA regressed out in trajectory models.

## 5. Promotion Proposals

| # | target skill | kind | summary | source |
|---|---|---|---|---|
| P1 | sc-preprocessing | requirement | For nuclei/snRNA input, count **exon+intron** (sum), never exon-only — introns carry most nuclear signal (parallels CellRanger `--include-introns`). | M2 |
| P2 | sc-annotation | method | Add **SVM 5-fold cross-validation F1 + permuted-label null** as a quantitative cluster-validity / over-clustering check (real≈0.99, permuted≈0.17). | M3 |
| P3 | sc-annotation | method | Add **Garnett marker-file classifier** (Monocle3) as an automated, literature-driven, portable annotation cross-check alongside SCimilarity/Azimuth/SingleR. | M4 |
| P4 | sc-integration | recipe | **Pan-compartment cross-tissue embedding**: extract one cell class from all organs, downsample ≤5000/type/organ, PCA on top-marker gene set, UMAP(cosine, n_neighbors=50)+Louvain(res=1e-4) — prevents the largest organ dominating. | M5 |
| P5 | sc-integration | recipe | **Cross-species bridge**: Seurat anchors on shared-ortholog HVGs (30 dims, 3000 HVGs) + k=3 NN label transfer + NNLS coeff<0.6 "no-1:1-match" gate. | M6 |
| P6 | scientific-plotting | recipe | **"Same embedding, recolored N ways"** panel set + **grey-background focal overlay** (context grey, focal set colored) + hand-drawn directional trajectory arrows. | F6, F5 |
| P7 | scientific-plotting | recipe | **Labeled per-organ UMAP small-multiple grid** with on-plot cluster text labels (no side legend) and a summary-stats text cell. | F1 |
| P8 | scientific-plotting | recipe | **Cross-validation confusion-matrix heatmap** (viridis, square tiles, diagonal-ordered) paired with a **specificity-F1 box plot vs permuted null**. | F2 |
| P9 | sc-conventions | requirement | **State exclusion thresholds in the figure/caption** (e.g. "samples ≤200 blood cells excluded") as a provenance rule for composition/proportion plots. | §4 |
| P10 | scientific-plotting | palette | Reinforce **Z-score-with-explicit-cap** ([0,3]) + **block-diagonal (argmax) ordering** as the default for marker expression heatmaps/dotplots. | F3, F4 |

`kind` ∈ { recipe, requirement, method, palette, fix }.

## 6. Directed distill tables

### 6.1 Paper info table

| title | authors | journal | year | doi | main findings summary |
|---|---|---|---:|---|---|
| A human cell atlas of fetal gene expression | Cao J, O'Day DR, Pliner HA, et al. (Trapnell & Shendure) | Science | 2020 | 10.1126/science.aba7721 | sci-RNA-seq3 atlas of ~4M cells across 15 human fetal organs; 77 main types / 657 subtypes; quantitative specificity framework (SVM-F1, Garnett); cross-organ analysis of blood/endothelial/epithelial; adrenal as minor erythropoietic site; trophoblast-like & hepatoblast-like cells in unexpected organs; straightforward human↔mouse (MOCA) integration showing conserved developmental programs. |

### 6.2 Marker/signature curation table

| marker_or_signature | marker_type | experiment | species | context | claim_use | evidence_source | citation |
|---|---|---|---|---|---|---|---|
| CD8B, CD5 | gene/protein | scRNA (sci-RNA-seq3) | human | fetal T cells (pan-organ) | Specific T-cell markers where canonical CD4/CD8A leak into mac/DC/NK | Fig 4D + main text | 10.1126/science.aba7721 |
| SORCS1, JMY (also RORC, KIT) | gene/protein | scRNA | human | fetal ILC3 (pan-organ) | ILC3-specific markers beyond RORC/KIT | Fig 4D + main text | 10.1126/science.aba7721 |
| OLR1, SIGLEC10, RP11-480C22.1 | gene/protein (incl. ncRNA) | scRNA | human | fetal microglia | Novel microglia markers alongside CLEC7A/TLR7/CCL3 | main text (microglia) | 10.1126/science.aba7721 |
| CD34, CD74, MPO | gene/protein | scRNA | human | HSPC (blood trajectory) | HSPC stage markers | Fig 5B | 10.1126/science.aba7721 |
| SLC16A9, FAM178B | gene/protein | scRNA | human | early erythroid progenitors (EEP) | Erythroid stage marker | Fig 5B | 10.1126/science.aba7721 |
| KIF18B, KIF15 | gene/protein | scRNA | human | committed erythroid progenitors (CEP) | Erythroid stage marker | Fig 5B | 10.1126/science.aba7721 |
| TMCC2, HBB | gene/protein | scRNA | human | erythroid terminal differentiation (ETD) | Erythroid stage marker | Fig 5B | 10.1126/science.aba7721 |
| KIT, LMO4, CPA3 | gene/protein | scRNA | human | Basophil/Mast | Lineage markers | Fig 5B | 10.1126/science.aba7721 |
| F13A1, RNASE1, COLEC12, LYVE1 | gene/protein | scRNA | human | perivascular macrophages | Macrophage subtype identity | Fig 5H + main text | 10.1126/science.aba7721 |
| TIMD4, CD5L, VCAM1 (MYO9A, NDST3) | gene/protein | scRNA | human | phagocytic macrophages | Macrophage subtype identity | Fig 5H + main text | 10.1126/science.aba7721 |
| SEMA6B, HLA-DPB1, HLA-DQA1, HLA-DPA1, AHR | gene/protein | scRNA | human | antigen-presenting macrophages (GI tract) | Macrophage subtype identity | Fig 5H + main text | 10.1126/science.aba7721 |
| TMEM119, CX3CR1 | gene/protein | scRNA | human | microglia (cerebrum-enriched) | Microglia subcluster marker | main text | 10.1126/science.aba7721 |
| SOX6, KLF1, E2F2 | signature (conserved TFs) | cross-species integration | human\|mouse | erythroblast identity | Human-mouse conserved erythroblast TFs | Fig 6E | 10.1126/science.aba7721 |
| PBX1, MEIS1, FHL2, MYLK, FLI1, MAFG | signature (conserved TFs) | cross-species integration | human\|mouse | megakaryoblast identity | Human-mouse conserved megakaryoblast TFs | Fig 6E | 10.1126/science.aba7721 |
| DAB2, TCF7L2, MAFB, NR1H3 | signature (conserved TFs) | cross-species integration | human\|mouse | macrophage identity | Human-mouse conserved macrophage TFs | Fig 6E | 10.1126/science.aba7721 |
| SALL1, MEF2A, NUAK1 | signature (conserved TFs) | cross-species integration | human\|mouse | microglia identity | Human-mouse conserved microglia TFs | Fig 6E | 10.1126/science.aba7721 |
| CSH1, CSH2 | gene/protein | scRNA + IF | human | trophoblast-like cells (lung, adrenal) | Circulating-trophoblast-like identity in unexpected organs | Fig 3A + fig S12 | 10.1126/science.aba7721 |
| AFP, ALB | gene/protein | scRNA + IF | human | hepatoblast-like cells (placenta, spleen) | Circulating-hepatoblast-like identity | Fig 3A,D + fig S12 | 10.1126/science.aba7721 |
| HMX1 (NKX5-3) | gene/protein | scRNA | human | chaffin cells / sympathetic neuron diversification | Chromaffin identity TF | main text (epithelial/neuroendocrine) | 10.1126/science.aba7721 |
| ASCL1, NKX2-1 | gene/protein | scRNA | human | pulmonary neuroendocrine cells (PNEC) | PNEC identity TFs | main text | 10.1126/science.aba7721 |
| HNF4A | gene/protein | scRNA | human | proximal tubule (kidney) | Proximal-tubule identity/required for formation | main text (epithelial) | 10.1126/science.aba7721 |
| MAFB, TCF21/POD1 | gene/protein | scRNA | human | podocytes (metanephric) | Podocyte identity TFs | main text | 10.1126/science.aba7721 |
| IGFBP1, DKK1 / PAEP, MECOM | gene/protein | scRNA + maternal genotype | human | maternal decidual stromal / endometrial epithelial (placenta) | Maternal-origin cell identity (XIST/TSIX high) | main text + fig S12B | 10.1126/science.aba7721 |
| STC2, TLX1 | gene/protein | scRNA | human | spleen mesenchymal precursor/stem-like | Annotation of initially unannotated type | main text + fig S13 | 10.1126/science.aba7721 |

## 7. Evidence level tags
- M1, M2 + code snippets: `repo-only` + `full-text` (methods) — pipeline code read directly at pinned SHA.
- M3–M7 (clustering, specificity, Garnett, cross-organ, cross-species, DDRTree): `methods-only` — described in STAR Methods; **no analysis/figure code shipped**.
- F1–F6: `captions-only` — every figure template is **reconstructed**, no plotting code in the repo.
- Markers (§6.2): `full-text` (main text + figure captions); some subtype panels are figure-label-derived.
