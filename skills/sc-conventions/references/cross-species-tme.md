# Cross-species TME analysis conventions
(source: https://doi.org/10.1038/s41590-026-02505-7)

Rules for any human–mouse tumor atlas. Plotting recipes live in `scientific-plotting/references/cross-species-tme.md`.

## Map mouse samples into a published human TME taxonomy
Do not stop at within-atlas proportion bars. Embed mouse genotypes / models in a published human
classification (ImmunoProfiler archetypes and/or Thorsson C1–C6) using a shared 8–12 feature
composition matrix, pooled z-score, and cosine similarity.
(source: https://doi.org/10.1038/s41590-026-02505-7)

Typical mice recapitulate immune-desert, macrophage-rich human TMEs (~17% of a pan-cancer
ImmunoProfiler cohort). State which human subset a GEMM actually models.

## Archetype-restricted sensitivity for frequency correlations
A conserved *abundance* is not a conserved *coupling*. If CD8–macrophage (or T–myeloid)
correlations appear in mouse, recompute them in human after restricting to desert / myeloid-rich
samples. Report both. Mouse-only couplings over-generalize to T-cell-dominant human TMEs.
(source: https://doi.org/10.1038/s41590-026-02505-7)

## Split T vs myeloid before NMF; match with Jaccard AND gene-weight scatter
Run cNMF / NMF separately on conventional T cells and non-granulocytic myeloid cells.
Match programs by Jaccard of top 20 (T) or top 50 (myeloid) genes with a Fisher exact test.
Confirm with a human-vs-mouse gene-weight scatter. Jaccard > 0.05 is not conservation.
Validate matched programs against external meta-programs (Gavish 2023; Courau ST2).
Drop contamination / doublet factors before interpreting movements.
(source: https://doi.org/10.1038/s41590-026-02505-7)

Granulocytes / PMN-MDSCs are outside the Courau myeloid GEP set — keep them as a separate claim.

## Score coordinated GEP movements, not isolated myeloid signatures
The conserved prognostic unit in Courau is T-cell cytotoxicity (T_3: PRF1, LAG3, GZMB, NKG7)
co-occurring with myeloid IFN response (My_2: IFIT2, IFIT3, ISG15, CXCL10). Stratify survival
as a 2×2 (T_3 Hi/Lo × My_2 Hi/Lo). Myeloid-high vs low is a different, often opposite, axis.
(source: https://doi.org/10.1038/s41590-026-02505-7)

Prefer scRNA-vs-scRNA over bulk-sorted-human vs scRNA-mouse; the latter underestimates conservation.

## Chemokine source cells, not just ligand–receptor pairs
CCC of conserved LR pairs does not test whether the ligand is made by the same compartment.
Z-score chemokine genes across compartments within species. Known failure mode: CXCL13 and
CXCL9/10 are T/myeloid-biased in human TMEs and stroma-biased in typical mouse models.
Conserved receptors to lean on: CXCR4, CXCR6, CCR4, CCR7.
(source: https://doi.org/10.1038/s41590-026-02505-7)
