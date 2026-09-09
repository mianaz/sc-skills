---
name: sc-target
description: Use when prioritising drug targets or disease-relevant cell types from scRNA-seq plus GWAS (scDRS, coloc, Open Targets).
---

# scRNA Disease & Drug Target Prioritisation

## When to use
- User wants to identify which cell types are most relevant to a GWAS-implicated disease
- User wants to score each cell for enrichment of GWAS heritability signal (scDRS)
- User wants to link GWAS loci to specific genes or cell types via colocalization
- User wants to aggregate multi-evidence target scores using Open Targets

## When NOT to use
- No annotated scRNA object → run `sc-annotation` first
- No GWAS summary statistics available → cannot run scDRS or colocalisation
- User is new to single-cell analysis → direct to `sc-preprocessing`, `sc-integration`, `sc-annotation` first; this is an advanced downstream step

## Workflow
1. **GWAS heritability enrichment per cell type:** Use LD score regression (LDSC) with cell-type-specific annotations derived from scRNA marker genes. Identifies which cell types are enriched for disease heritability.
2. **Per-cell disease scoring (scDRS):** Given GWAS summary stats and a scRNA object, scDRS scores each cell for enrichment of the top GWAS genes. Produces a cell-level disease score; visualise on UMAP.
3. **Colocalization (coloc/coloc2):** For GWAS loci overlapping eQTLs (from GTEx or cell-type-specific eQTL studies), test whether the GWAS signal and eQTL signal share a causal variant. Colocalising loci nominate the causal gene at each locus.
4. **Open Targets aggregation:** query the public Open Targets GraphQL API for target–disease association scores, genetic evidence, and tractability; cross-reference with your colocalization and scDRS results to rank targets. Self-contained call (no external skill needed):
   ```python
   import requests
   q = """query($ensg:String!,$efo:String!){
     disease(efoId:$efo){ associatedTargets(BFilter:$ensg){ rows{ target{approvedSymbol} score } } } }"""
   r = requests.post("https://api.platform.opentargets.org/api/v4/graphql",
                     json={"query": q, "variables": {"ensg": "ENSG00000141510", "efo": "EFO_0000305"}})
   r.json()   # association score + evidence; see platform-docs.opentargets.org for the schema
   ```
   (The `database-lookup` skill is an optional convenience wrapper for this and GTEx eQTLs — not required.)
5. **Cancer dependency validation (optional):** cross-reference top targets against DepMap CRISPR (Chronos) dependency scores to check essentiality in relevant lineages — download the gene-effect matrix from depmap.org (or the `depmap` skill, optional) and look up your targets' scores.

## Key caveats
- scDRS results depend heavily on GWAS gene set quality (top N genes); always test sensitivity to N
- Colocalization requires fine-mapped or well-tagged loci; poorly tagged GWAS signals produce unreliable coloc posteriors
- Open Targets scores aggregate evidence but are not a substitute for experimental validation
- This workflow produces a prioritised candidate list, not validated targets

## Routes to
- `sc-grn` — TF regulatory context of prioritised target genes (internal suite)
- `database-lookup` *(optional, external)* — convenience wrapper for Open Targets / GTEx eQTL lookups
- `depmap` *(optional, external)* — CRISPR dependency lookups for top targets
- `population-genetics` *(optional, external)* — GWAS-to-function steps (TWAS, coloc, MR)
