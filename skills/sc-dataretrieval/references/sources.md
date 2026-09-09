# Finding & downloading public scRNA-seq data

Goal: get the deposited count files (preferred) or fastq (last resort) for an accession, plus a copy of the authors' processing description. Record download date and any checksums in methods.md.

## GEO (GSE / GSM)
GEO holds the **supplementary files** authors uploaded — this is where deposited count matrices live.

- Series landing page: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXX`
  - "Supplementary file" section (series-level) often has a `GSEXXXXXX_RAW.tar` bundling all samples.
  - Each sample (GSM) page lists its own supplementary files — typically the 10x triplet (`barcodes.tsv.gz`, `features.tsv.gz`/`genes.tsv.gz`, `matrix.mtx.gz`), an h5, a per-sample count csv, or a processed object.
- Bulk download of all supplementary files for a series:
  ```bash
  # FTP mirror path: last 3 digits of the accession become "nnn"
  # GSE164nnn for GSE164522, etc.
  wget -r -np -nH --cut-dirs=4 -A '*' \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164522/suppl/"
  ```
- Programmatic listing of files (so you can see tiers before downloading):
  ```bash
  # GEOquery (R) — metadata + supplementary file URLs, no bulk download
  Rscript -e 'library(GEOquery); g <- getGEO("GSE164522", GSEMatrix=FALSE); head(Meta(g))'
  Rscript -e 'GEOquery::getGEOSuppFiles("GSE164522", fetch_files=FALSE)'  # lists supp file URLs
  ```
- Inspect a `_RAW.tar` without extracting everything: `tar -tvf GSEXXXXXX_RAW.tar`.
- The GSM/GSE pages' **"Data processing"** and **"Extraction/Library"** fields are the documented preprocessing — see references/provenance.md.

## ArrayExpress / BioStudies (E-MTAB-…, E-GEOD-…)
ArrayExpress data is served via BioStudies.

- Landing page: `https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-XXXX`
- File list JSON (lists every processed/raw file with download URLs):
  ```bash
  curl -s "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-XXXX" | less   # find "files" section
  ```
- Processed files (count matrices, objects) are under the study's `Files` → download directly via the listed FTP/HTTP URLs.
- The **IDF/SDRF** sample-and-data-relationship files document protocols and processing; grab them for provenance.

## SRA / ENA (fastq — tier 4 only)
Only when no usable counts are deposited. fastq → CellRanger → raw matrix → CellBender.

- Find the run accessions (SRR…/ERR…) from the GEO/ArrayExpress sample pages or:
  ```bash
  # ENA: get a TSV of fastq URLs + md5 for a study/sample
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRPXXXXXX&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
  ```
- Download fastq from the ENA `fastq_ftp` URLs (resumable, no SRA toolkit needed); verify `fastq_md5`.
- Or use `prefetch`/`fasterq-dump` (SRA Toolkit) if you must pull from NCBI SRA.
- 10x fastqs must keep CellRanger naming (`SAMPLE_S1_L001_R1_001.fastq.gz` etc.) or be renamed before `cellranger count`.

## Sanity checks before loading
- Confirm gene IDs vs symbols, and organism, match what you expect.
- Confirm cell count per sample is plausible (raw ≫ filtered).
- Note per-sample which file you took and its tier — this drives load-by-format.md.
