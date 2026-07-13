# scVI / scANVI (Python, adaptable path)

Export Seurat → h5ad with counts in a layer (see sc-conventions r-python-bridge).
```python
import scvi, scanpy as sc
adata = sc.read_h5ad("objects/seu.h5ad")
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")
model = scvi.model.SCVI(adata, n_latent=30)
model.train(max_epochs=400)                      # adjust by dataset size + convergence
adata.obsm["X_scVI"] = model.get_latent_representation()

# scANVI: semi-supervised, needs seed labels in adata.obs["seed_labels"]
scanvi = scvi.model.SCANVI.from_scvi_model(
    model, adata=adata, labels_key="seed_labels", unlabeled_category="Unknown")
scanvi.train(max_epochs=200)
adata.obsm["X_scANVI"] = scanvi.get_latent_representation()
adata.obs["scanvi_pred"] = scanvi.predict()

import pandas as pd
pd.DataFrame(adata.obsm["X_scVI"],   index=adata.obs_names).to_csv("objects/X_scVI.csv")
pd.DataFrame(adata.obsm["X_scANVI"], index=adata.obs_names).to_csv("objects/X_scANVI.csv")
adata.obs[["scanvi_pred"]].to_csv("objects/scanvi_pred.csv")
```
Carry back into Seurat via CreateDimReducObject (see r-python-bridge): the scVI latent
into `seu[["scvi"]]` (key `scVI_`) and, if you used the scANVI path, the scANVI latent
into `seu[["scanvi"]]` (key `scANVI_`). Cluster/UMAP on the chosen reduction
(`"scvi"` or `"scanvi"`), and carry `scanvi_pred` into `meta.data`.
Log n_latent, batch_key, epochs, rationale.
