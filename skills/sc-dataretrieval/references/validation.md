# Validating an accession before you trust it

**The deposited label lies more often than you'd think.** In a 298-accession liver curation pass, validating every entry against GSM sample-level metadata + the paper found that **~55% of "conflict" rows had a wrong label**: wrong tissue, wrong assay (bulk mislabeled single-cell), wrong species, or sample-count ≠ individual-count. Series-level `title` / `gdsType` and any auto-curated `data_kind` are **hypotheses, not facts** — confirm them from the GSM level before downloading or analyzing.

Do this BEFORE format-tiering (sources.md) — there's no point picking a file tier for a series that turns out to be kidney ATAC-seq.

## The five things to confirm (all from GSM-level data, not the series blurb)

| Field | How to confirm | Failure mode it catches |
|---|---|---|
| **Species** | Tally `taxon` across GSMs, not the series taxon | "human" series that is human+mouse (humanized mice, xenografts, species-mixing controls) |
| **Assay / data_kind** | `gdsType` + GSM **protocol/files** (not the title) | `genome binding/occupancy` = ATAC/ChIP/CUT&RUN; titles `RNA-Seq`/`WES`/`minibulk`/`lpWGS`/`eCLIP`/`MPRA` = **not single-cell**; "Multiome/Spatial/snRNA" labels are frequently bulk subseries; **"spatial" is often an in-silico analysis of scRNA, not a spatial assay → see "The spatial trap" below** |
| **Tissue** | GSM titles + paper | "liver" series that is PBMC/blood, adipose, organoid, or a cell line (HepG2/Huh7/Hep3B…); disease-relevant ≠ tissue-of-origin |
| **n_samples vs n_individuals** | Count GSMs; parse donor IDs from titles | replicates, timepoints, sorted fractions, plate wells, and pooled donors all inflate GSM count above true biological-individual count |
| **PMID ↔ dataset match** | Compare paper abstract subject to GEO subject | sheets/links sometimes cite the **wrong paper** for an accession |

## The "spatial" trap: assay vs analysis

**A title containing "spatial" does NOT mean a spatial assay.** A large share of liver
"spatial" studies are scRNA/snRNA-seq that *computationally reconstruct* lobule zonation
("Spatial Reconstruction of…", "spatially-sorted hepatocytes", "spatial transcriptional
activity", zonation mapping). The deposited data is a gene×cell matrix — single-cell —
with no spatial assay at all. Confirmed case: GSE200771 ("Spatial Reconstruction of the
Early Hepatic Transcriptomic Landscape…") is 10x Chromium 3′ scRNA-seq; the word "spatial"
is an in-silico analysis, not Visium. Classifying it as `data_kind=Spatial` is wrong, and
calling it "verified" is worse — it launders a title keyword into a fact.

**Decide assay from mechanics + files, never from the title:**

| Signal | Real spatial **assay** | Spatial **analysis** of scRNA (→ label scRNA-seq) |
|---|---|---|
| Platform / protocol | Visium, Xenium, MERFISH, Slide-seq, GeoMx, CosMx, seqFISH, Stereo-seq, spaceranger | "single-cell suspension", 10x Chromium 3′/5′, CellRanger, 16bp 10x barcode + UMI |
| Deposited files | `tissue_positions*`, `scalefactors_json`, fiducial/H&E `.tif`, `spatial/` dir | gene×cell matrix (`matrix.mtx` + `barcodes`/`features`), per-cell h5/h5ad |
| Title words that DON'T prove spatial | — | "spatial reconstruction/mapping", "zonation", "spatially-sorted", "spatial activity/redistribution" |

Watch two failure modes when checking: (1) **sample-cap** — if you only read the first
N GSMs, a multi-modal study's spatial samples may sit *after* the cap; scan all GSM titles
for the spatial-platform keywords before concluding "no spatial". (2) **over-promotion** —
a true Visium study often says "single-cell suspension" in its dissociation step; that lone
phrase is NOT enough to add "single-cell" to a correct `Spatial` label. Promote to
`Single-cell; Spatial` only when a distinct single-cell **assay** (its own 10x GEX GSMs) is
present, not merely sc-sounding protocol words.

## The superseries trap
An accession may resolve to **one subseries** of a superseries while the data you want lives in a sibling. Example seen: a series record resolved to the `[RNA-seq]` cell-line subseries (52 bulk samples) while the patient scRNA (n=15) and spatial (n=22) lived in two *different* GSE numbers. Always check whether the GSM set you got matches the assay/n you expect; if not, look for sibling subseries.

## How to pull GSM-level metadata (no API key needed; ~3 req/s anon, 10 with key)

NCBI E-utilities on `db=gds`:
```bash
EUTILS=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
# 1. series UID
curl -s "$EUTILS/esearch.fcgi?db=gds&term=GSE242889%5BACCN%5D+AND+gse%5BETYP%5D&retmax=1"
# 2. series summary (title, taxon, gdsType, n_samples, and nested Samples list of GSMs)
curl -s "$EUTILS/esummary.fcgi?db=gds&id=<UID>"
# 3. per-GSM taxon: esearch each GSM accession, esummary the UID, read its <Item Name="taxon">
```
Or via GEOquery in R (same idea, returns the per-sample `taxon`/`characteristics_ch1`/`title`):
```r
g <- GEOquery::getGEO("GSE242889", GSEMatrix=FALSE)
sapply(GEOquery::GSMList(g), function(s) GEOquery::Meta(s)$organism_ch1)   # species per GSM
sapply(GEOquery::GSMList(g), function(s) GEOquery::Meta(s)$title)          # titles → donors/conditions
```

## Cross-check the publication
- **EuropePMC** for open-access Methods + Data-availability:
  ```bash
  EPMC=https://www.ebi.ac.uk/europepmc/webservices/rest
  curl -s "$EPMC/search?query=EXT_ID:37972953+AND+SRC:MED&resultType=core&format=json"  # → pmcid, isOpenAccess
  curl -s "$EPMC/PMC1234567/fullTextXML"                                                  # OA only
  ```
- The abstract usually states the **true individual count** ("ten specimens from five patients" → n_individuals=5, not 10) and the disease/tissue.
- **Paywalled = you cannot fully verify.** Mark the field `unverified`; never invent a value to fill a cell.

A ready-made helper that does all of the above for one accession lives at `tmp/liver_db/fetch_meta.py` (esummary + per-GSM taxon tally + GSM titles + PubMed abstract + EuropePMC OA methods; self-throttles).

## Record the verdict
For each accession capture: `validation_status` ∈ {`verified`, `conflict`, `unverified`}, `geo_verified_fields` (what you actually confirmed from GEO), and a one-line note stating *sheet-value vs found-value* for any conflict. This is what lets a human (or the next agent run) trust the row without re-deriving it.

## What stays a human decision
Validation produces facts; **inclusion is policy.** "Is blood from a liver-disease cohort in scope?" / "Is a HepG2 cell-line dataset a 'liver' dataset?" are judgment calls the agent should *flag and tier*, not silently drop. Keep clearly-out (wrong tissue, bulk) separable from borderline (mets, hepatic cell line, tissue-is-blood) so a person sets the boundary.
