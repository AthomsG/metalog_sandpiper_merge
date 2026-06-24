# Metalog × Sandpiper merge (human / ocean)

Find the sequencing runs that exist in **both** metalog and sandpiper for the
**human** and **ocean** environments, and count species occurrences in them.

## Data

Downloaded by `01_fetch_data.py` into `fetched_data/` (from metalog.embl.de):

- `human_extended_wide_latest.tsv.gz` — all human-environment samples (metadata, wide).
- `ocean_extended_wide_latest.tsv.gz` — all ocean-environment samples.
- `sequencing_db_mapping_latest.tsv.gz` — maps metalog samples to runs/experiments/sample accessions.

Sandpiper (also fetched by `01_fetch_data.py`, into `/lisc/data/work/.../sandpiper/`):
the **latest** version is resolved from the Zenodo concept record and the GTDB profile
file (e.g. `sandpiper2.0.1.gtdb.csv.gz`, ~3.3 GB) is downloaded — taxonomic profiles,
one row per (run, taxon), keyed by **run accession**. Because the version is in the
filename, the big file is only re-downloaded when a **newer version** is published
(otherwise the existing local copy is reused after a size check; fresh downloads are
MD5-verified). The resolved path is written to `fetched_data/sandpiper_latest.txt`, which
`02_merge_and_count.py` reads. (Earlier we used the GlobDB profile `sandpiper1.1.0.globdb.csv.gz`;
switching taxonomy is a one-line `SANDPIPER_PATTERN` change.)

## Samples vs. runs (important)

- Metalog organizes things by **sample** (its own `sample_alias`, e.g. `Cait_2019_infant.sample_AC114`).
- Sandpiper is keyed by **run accession** (e.g. `SRR11852051`, `DRR000713`).
- These are not the same ID, so we use the mapping file to join them:

  ```
  wide.sample_alias  →  mapping (kind=="run")  →  mapping.external_id (run accession)  →  sandpiper.sample
  ```

- **One sample can have multiple runs** (up to ~35 here). A sample is one biological
  specimen; a run is one sequencing run of it (re-sequencing, lanes, replicates...).
  We keep **all** runs of every sample. So the run counts are larger than the sample
  counts, and a species detected in several runs of the same sample is counted once
  per run.

## Run

```
python3 01_fetch_data.py       # metalog files (always fresh) + latest sandpiper (if newer)
python3 02_merge_and_count.py  # single pass over sandpiper, handles both environments
```

## Output

`output/human/` and `output/ocean/`, each containing:

- `shared_runs.txt` — run accessions present in both that environment's metalog set and
  sandpiper (one per line).
- `species_occurrences.csv` — `species,occurrences`; `occurrences` = number of sandpiper
  rows among those shared runs in which the `s__` species appears (≈ number of shared
  runs detecting it).

## Coverage threshold

A species is only counted in a row whose **coverage ≥ `MIN_COVERAGE`** (default `10`),
set at the top of `02_merge_and_count.py`. This filters `species_occurrences.csv` only;
`shared_runs.txt` lists every run shared between metalog and sandpiper regardless of
coverage.
