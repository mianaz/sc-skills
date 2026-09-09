# Topology-faithful embedding (PAGA-initialized force-directed graph)

When you need a 2D layout that preserves the **branching topology** of differentiation (UMAP/t-SNE tear or merge connections), use a force-directed graph (ForceAtlas2) **seeded by PAGA** — not UMAP alone. PAGA gives coarse cluster connectivity; `init_pos="paga"` plants the single-cell FA2 layout on it so global branch structure survives while fine structure spreads out. This is the classic haematopoiesis layout.

(source: sc-paper-distill/papers/2019-FCA-liver.md M1, M2, F4 — Popescu et al. Nature 2019; AGA in that paper is the precursor to PAGA.)

```python
import scanpy as sc
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.leiden(adata)
sc.tl.paga(adata, groups="leiden")          # coarse connectivity (AGA's successor)
sc.pl.paga(adata, plot=False)               # must run before draw_graph init
sc.tl.draw_graph(adata, init_pos="paga", layout="fa")   # ForceAtlas2; needs the fa2 package
sc.pl.draw_graph(adata, color="leiden")
# orient branches: sc.tl.dpt(adata) with an HSC/progenitor-marker root cell
```

`layout="fa"` needs `fa2`; fall back to `layout="fr"` (Fruchterman–Reingold) if unavailable. The `init_pos="paga"` chaining is the load-bearing step — `draw_graph` alone does not preserve topology. Diffusion-map pseudotime here requires a **single true lineage** (remove doublets/outliers first) or it produces nonsense.
