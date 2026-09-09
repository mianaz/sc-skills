# scGPT (foundation-model embeddings)

scGPT is a pretrained transformer for single-cell data. In this suite it is a **Python/GPU detour**: embed cells → carry the embedding back into the Seurat object (sc-conventions r-python-bridge) → cluster/annotate as usual. Use it as an alternative embedding (like scVI) and as one candidate-annotation signal — never as the final label (markers decide, per sc-annotation).

## When to use
- Foundation-model cell embeddings for clustering/integration (alternative to scVI/Harmony)
- Embedding-space annotation: kNN/label-transfer from a labelled reference embedded the same way
- Gene-level representations for perturbation / GRN context

## When not to use
- Probabilistic batch-corrected latent + Bayesian DE → scVI/scANVI (`references/scvi-scanvi.md`)
- Final identity assignment → never from scGPT alone; verify with markers (sc-annotation)
- No GPU (≥16 GB VRAM) → Harmony or scVI instead

## Requirements
`scgpt` (Python ≥3.10, CUDA 12.1+), GPU ≥16 GB VRAM (24 GB+ for large data). The checkpoint is a **raw directory** (`args.json`, `best_model.pt`, `vocab.json`), ~200 MB for the released human model — download once and reference the path. Run as a batch/GPU job, not on a login node (scientific-reproducibility).

## Embed → carry back to Seurat
```python
import anndata as ad, pandas as pd
from scgpt.tasks import embed_data
adata = ad.read_h5ad("objects/seu.h5ad")             # exported from Seurat (sc-conventions bridge)
emb = embed_data(adata, model_dir="/path/scgpt-human",
                 gene_col="feature_name",            # adata.var column holding gene symbols
                 use_fast_transformer=False)         # True only if flash_attn imports cleanly
pd.DataFrame(emb.obsm["X_scGPT"], index=emb.obs_names).to_csv("objects/scgpt_emb.csv")  # n×512
```
```r
z <- as.matrix(read.csv("objects/scgpt_emb.csv", row.names = 1)); z <- z[colnames(seu), ]
colnames(z) <- paste0("scGPT_", seq_len(ncol(z)))
seu[["scgpt"]] <- Seurat::CreateDimReducObject(z, key = "scGPT_", assay = DefaultAssay(seu))
# cluster/UMAP on reduction = "scgpt" as usual
```
Genes absent from the scGPT vocabulary are dropped — point `gene_col` at the right symbol column and standardize symbols first (sc-preprocessing/gene-symbols).

## Annotation from the embedding
The released model ships embeddings, not an annotation head. For labels: embed a labelled reference with the same `embed_data`, then kNN/label-transfer in scGPT space; or cluster the embedding and label clusters by markers. Either way it is ONE candidate signal → confirm with sc-annotation's mandatory marker dotplot.
