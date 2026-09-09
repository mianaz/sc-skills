# Marker-panel discriminatory power (QC, not production labelling)

Quantify whether a chosen gene panel actually *separates* your annotated cell types — e.g. before
committing to a flow-cytometry panel. Frame it as supervised classification on panel genes only,
and benchmark against sensible ceiling/floor. Diagnostic metric, not a way to assign labels.
(source: sc-paper-distill/papers/2019-FCA-liver.md M3 — Popescu et al. Nature 2019, "gene discriminatory power", flow-gating analogy.)

## Recipe

```python
import scanpy as sc, numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score

panel = [...]                                   # your candidate marker genes
X = adata[:, panel].X.toarray()
y = adata.obs["cell_type"].values
groups = adata.obs["donor"].values              # CV by donor → no batch/identity leakage

clf = RandomForestClassifier(n_estimators=500, class_weight="balanced", n_jobs=-1)
pred = cross_val_predict(clf, X, y, groups=groups, cv=GroupKFold(5))
print(classification_report(y, pred))           # per-class precision/recall/F1 (not just accuracy)
cm = confusion_matrix(y, pred, labels=np.unique(y))   # off-diagonal = confusable type pairs
```

## Make the number meaningful

- **Ceiling & floor:** rerun on the full HVG set (ceiling) and on a random/size-matched gene set
  or permuted labels (floor). "Sufficient panel" = close to ceiling, far above floor.
- **Per-class, not global:** report macro-F1 / balanced accuracy and one-vs-rest ROC- or PR-AUC
  (PR-AUC for rare populations). The **confusion matrix is the key output** — it names the cell-type
  pairs the panel cannot separate (the flow-panel failure mode).
- **CV by donor/sample** (GroupKFold), never plain k-fold, or you leak identity across folds.
- **Feature importance** (RF importances / SHAP) tells you which markers are load-bearing and which
  are redundant and can be dropped.

## Flow-panel refinements

- Restrict the panel to **surface / antibody-available** markers before training.
- Flow gates on a few channels, so also check **pairwise separability** (2D density / shallow trees
  on marker pairs) — confirm each target is isolable by an achievable gating sequence, not only by a
  high-dimensional boundary.

Validated cell-type labels (the `y`) must come from the annotation flow first; this only tests a panel.
