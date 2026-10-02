# Cross-species integration & label transfer

Bridging datasets across species (human ↔ mouse, and further) is mostly a **shared-ortholog + anchor**
problem plus a **principled "no match" gate** so you don't force a 1:1 mapping that doesn't exist.
(source: https://doi.org/10.1126/science.aba7721)

## Recipe

1. **Map to one-to-one orthologs** and rename genes to a common symbol space (e.g. biomaRt / gprofiler2
   ortholog table; keep only 1:1 orthologs to avoid many-to-one collisions).
2. **Co-embed with Seurat v3 anchors** on shared-ortholog HVGs.

```r
shared_hvg <- intersect(VariableFeatures(human), VariableFeatures(mouse))   # in common symbols
anchors <- Seurat::FindIntegrationAnchors(list(human, mouse), anchor.features = shared_hvg,
                                          dims = 1:30)                        # Cao: 30 dims, ~3000 HVGs
comb    <- Seurat::IntegrateData(anchors, dims = 1:30)
comb    <- Seurat::ScaleData(comb) |> Seurat::RunPCA() |> Seurat::RunUMAP(dims = 1:30)
```

3. **Transfer labels by k-NN in the joint space** with a *small* k so rare types stay resolvable
   (Cao used **k = 3**): each query cell takes the majority label of its k nearest reference-species cells.
4. **Gate no-1:1-match with NNLS.** Regress each query type's mean profile on the reference type
   means; a best-fit **coefficient < ~0.6 means "no confident match"** — report it as unmatched rather
   than assigning the nearest label. (Cao flagged mouse placenta/skin/gonads this way.)

## Why these choices

- **Shared-ortholog HVGs**, not all genes: species-divergent genes inject noise; conserved variable
  genes carry the cell-identity signal that actually aligns.
- **Small k for transfer:** large k averages rare populations into their abundant neighbours.
- **The NNLS gate is the point of the method** — cross-species mapping's main failure mode is
  confidently mislabelling a population that simply has no counterpart in the other species. Always
  report the no-match set explicitly (a `log()`-style note), never silently drop or force it.
- Expect **stage/age offset**: developmental datasets project onto the *nearest* time point of the
  other species (human fetal → late mouse embryonic in Cao) — interpret matches as cell-type, not age, equivalence.

Log ortholog source, HVG count, dims, k, and the NNLS threshold in methods.md.
