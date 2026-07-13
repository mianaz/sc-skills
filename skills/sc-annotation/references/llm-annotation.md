# LLM-based candidate annotation (GPT-4o / Claude via API)

Feed each cluster's top +/- markers (and optionally enriched gene sets) to an LLM; get a
cell type + confidence per cluster. **Candidates only** — confirm with the marker dotplot gate.
Works against any OpenAI-compatible chat endpoint (OpenAI, OpenRouter → Claude, etc.).

## Inputs
- `de_markers`: a FindAllMarkers data.frame (cluster, gene, avg_log2FC, p_val_adj).
- Optional `gsea_data`: per-cluster enriched pathways (named list keyed by cluster) — see the
  GSEA helper in sc-annotation (or decoupleR) to produce it.
- API key from an env var (`OPENAI_API_KEY` / `OPENROUTER_API_KEY`); never hardcode.

## Pattern (per cluster)
1. Filter markers: `p_val_adj < 0.05`, drop noise genes (e.g. mouse `Rik$|^Gm`).
2. Take top-N positive and top-N negative genes by avg_log2FC.
3. Build one prompt listing, per cluster, positive markers / negative markers / enriched pathways.
4. System prompt fixes the output format: `Cluster X: Cell_Type (0.xx)` with a confidence rubric;
   instruct "be specific with subtypes; never predict 'unknown'".
5. POST with retry/back-off; parse lines with a regex into `cluster | predicted_label | confidence`.

## Minimal call shape
```r
# call_llm_api(): POST {model, messages:[system,user]} with Authorization: Bearer <key>,
#   retry 3x with exponential back-off, extract choices[[1]]$message$content robustly.
res <- run_aiCelltype(seu, de_markers = markers, species = "mouse", tissue = "prostate",
                      gsea_data = gsea, top_n = 20,
                      base_url = "https://api.openai.com/v1/chat/completions",
                      model_name = "gpt-4o", api_key = Sys.getenv("OPENAI_API_KEY"))
# Claude via OpenRouter: base_url="https://openrouter.ai/api/v1/chat/completions",
#   model_name="anthropic/claude-sonnet-4.5", api_key=Sys.getenv("OPENROUTER_API_KEY")

labels <- setNames(res$predicted_label, as.character(res$cluster))
seu$celltype_llm <- labels[as.character(seu$seurat_clusters)]
```

## Output-parsing regex (robust to provider formatting)
```r
m <- regexec("^Cluster (\\d+): (.+?)\\s*\\((0?\\.[0-9]+|1\\.0?)\\)\\s*$", line, perl = TRUE)
```

Run two independent models (e.g. GPT-4o + Claude) and compare; agreement raises confidence.
All LLM labels are CANDIDATES — the canonical-marker dotplot (marker-verification.md) is still
the mandatory gate. Log model name, prompt design, and per-cluster confidence in methods.md.
A full reference implementation of `call_llm_api()` + `run_aiCelltype()` + a `run_gsea()` helper
(fgsea/clusterProfiler over msigdbr collections) is the pattern this is distilled from.
