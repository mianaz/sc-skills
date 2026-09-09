---
name: sc-trajectory
description: Use when inferring differentiation trajectories, pseudotime, RNA velocity, or cell fate from scRNA-seq.
---

# Single-Cell Trajectory & Pseudotime

## When to use
- Order cells along a differentiation path without time labels → pseudotime
- Spliced/unspliced counts and directional arrows on UMAP → RNA velocity (scVelo)
- Fate probabilities or terminal-state absorption → CellRank
- Graph-level connectivity before committing to a trajectory → PAGA

## When NOT to use
- Discrete unconnected compartments (e.g. lymphocytes vs epithelial) → skip
- No continuous biological process (e.g. cross-tissue atlas) → skip
- Raw counts not yet cleaned → `sc-preprocessing` first

## Method decision tree
1. Do cells form a continuous gradient on UMAP? No → stop.
2. Spliced/unspliced counts (STARsolo, alevin, kallisto|bustools with lamanno)? Yes → scVelo.
3. Directional fate probabilities from velocity? Yes → CellRank on top of scVelo.
4. No spliced/unspliced, or velocity unreliable? → PAGA or diffusion pseudotime.
5. Graph-level connectivity first? → PAGA, then refine.

## Layouts and R-native options
- Branching topology that UMAP tears → [references/paga-force-layout.md](references/paga-force-layout.md)
- Destiny / principal curve / slingshot in R → [references/r-pseudotime.md](references/r-pseudotime.md)

## Key caveats
- Pseudotime is a relative manifold ordering, not clock time.
- RNA velocity assumes constant splicing/degradation rates — sanity-check arrows against known biology.
- Wrong root inverts pseudotime; validate with progenitor markers.
- CellRank fates are only as good as the upstream velocity.

## Routes to
- `scvelo` — RNA velocity and latent time
- `scanpy` — PAGA and diffusion pseudotime
