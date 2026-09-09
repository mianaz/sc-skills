# Cluster specificity (is this a real cell type, or over-clustering?)

Turn "is this cluster a bona-fide, separable cell type?" into a **number with a null**, instead of
eyeballing marker separation. Use it to decide keep-vs-merge and to defend annotation granularity.
Diagnostic on the *clustering*, not a way to assign labels.
(source: sc-paper-distill/papers/2020-descartes-fetal-atlas.md M3 — Cao et al. Science 2020, SVM cross-validation specificity score.)

This complements `panel-discriminatory-power.md`: that tests whether a **small marker panel** separates
types; this tests whether a cluster/annotation is separable on the **whole transcriptome**.

## Recipe

```python
import scanpy as sc, numpy as np
from sklearn.svm import LinearSVC
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import f1_score

# subsample to balance types (Cao: ~2000 cells/type) so F1 isn't driven by the biggest cluster
def specificity_f1(adata, label_key, n_per=2000, n_splits=5, seed=0):
    idx = (adata.obs.groupby(label_key)
           .apply(lambda d: d.sample(min(len(d), n_per), random_state=seed)).index
           .get_level_values(-1))
    sub = adata[idx]
    X, y = sub.X, sub.obs[label_key].values
    pred = cross_val_predict(LinearSVC(), X, y, cv=StratifiedKFold(n_splits, shuffle=True, random_state=seed))
    return f1_score(y, pred, average=None, labels=np.unique(y)), np.unique(y)

f1, types = specificity_f1(adata, "cell_type")            # per-type CV F1 = specificity score
rng = np.random.default_rng(0)                            # permuted-label NULL (floor)
adata.obs["perm"] = rng.permutation(adata.obs["cell_type"].values)
f1_null, _ = specificity_f1(adata, "perm")
```

## Read it

- **Real main types score high** (Cao: median F1 ≈ 0.99), **finer subtypes lower but still well above the
  permuted null** (subtypes ≈ 0.77 vs permuted ≈ 0.17). A cluster whose F1 collapses toward the permuted
  floor is **over-clustering** — merge it into its nearest neighbour.
- **Always report the permuted-label null**, not just the raw F1: the null is the "what would random
  labels score" floor that makes a given F1 interpretable (and it absorbs class-imbalance/depth quirks).
- **Per-type, not global.** The point is *which* clusters are unstable; a mean F1 hides them.
- **Whole transcriptome, linear SVM, subsample per type.** A linear kernel keeps it a separability test
  (not a capacity test); per-type subsampling stops the largest cluster from dominating.
- Same procedure adjudicates **merge decisions** (score before vs after merging candidate pairs) and
  **subtype validity** within a compartment.

Feeds naturally into a confusion-matrix + F1-boxplot figure — see
`scientific-plotting/references/heatmaps.md` (cross-validation confusion matrix).
