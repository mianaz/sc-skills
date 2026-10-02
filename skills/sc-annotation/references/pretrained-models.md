# Pretrained annotation models (candidates only)

## Azimuth (R)
```r
seu <- Azimuth::RunAzimuth(seu, reference = "pbmcref")   # pick the matching reference
# predictions land in seu@meta.data: predicted.celltype.l1 (coarse),
# predicted.celltype.l2 (fine) + predicted.celltype.l2.score
# other references: lungref, kidneyref, bonemarrowref, heartref, pancreasref,
# tonsilref, adiposeref, fetusref, humancortexref ...
```

## SCimilarity (Python) — carry predictions back
Real API (module `scimilarity.cell_annotation`, class `CellAnnotation`):
```python
from scimilarity import CellAnnotation
from scimilarity.utils import align_dataset
import scanpy as sc

ca = CellAnnotation(model_path="/path/to/scimilarity_model")  # downloaded model dir
adata = sc.read_h5ad("counts.h5ad")                            # raw counts
adata = ca.annotate_dataset(adata)        # adds obs["celltype_hint"] + obsm["X_scimilarity"]

# export predicted labels keyed by barcode, read back into meta.data (sc-conventions bridge)
adata.obs[["celltype_hint"]].to_csv("scimilarity_pred.csv")
```
Lower level: `embeddings = ca.get_embeddings(align_dataset(adata, ca.gene_order).X)`
then `predictions, nn_idxs, nn_dists, stats = ca.get_predictions_knn(embeddings)`
(`stats` has per-cell confidence: hits, min_dist, vs2nd, vsAll). Read the CSV back into
`seu@meta.data` keyed by barcode (sc-conventions bridge).

## Garnett (R) — marker-file classifier (trainable, portable)
Unlike embedding models (SCimilarity/Azimuth) that match to a fixed reference, Garnett trains a
classifier from a **hand-written marker file** — so it is literature-driven, auditable, and
independent of your own clustering (a genuine cross-check, not a circular one).
(source: https://doi.org/10.1126/science.aba7721)

```r
library(garnett); library(org.Hs.eg.db)   # Garnett-for-Monocle3

# marker file: one ">CellType" block per type, with expressed: GENE1, GENE2 (from the literature)
classifier <- train_cell_classifier(
  cds = cds, marker_file = "markers.txt", db = org.Hs.eg.db,
  cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")
cds <- classify_cells(cds, classifier, db = org.Hs.eg.db, cluster_extend = TRUE,
                      cds_gene_id_type = "SYMBOL")   # -> colData(cds)$cell_type / cluster_ext_type
```

- **Why it earns its place:** because the marker file is written before seeing clusters, Garnett vs
  manual concordance is a real validity signal; a trained classifier is also **portable** — apply the
  same model to new datasets (Cao et al. validated a fetal-pancreas model on adult-pancreas scRNA-seq).
- `cluster_extend = TRUE` propagates confident calls to their cluster neighbours; report per-type
  concordance with your manual labels, and treat disagreements as clusters to re-examine, not to overwrite.
- Marker file quality is everything — draw genes from `sc-annotation/references/markers.md` (cited).

All predictions are CANDIDATES — confirm with markers before trusting.
