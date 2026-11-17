# DECISIONS

## 2026-11-17 18:30 UTC — OBIS source: AWS S3 GeoParquet primary, API fallback

**Decision:** Configure OBIS downloads to use AWS S3 GeoParquet as the primary source (via DuckDB), with automatic fallback to the API.

**Rationale:**
- Both `mammals_marine` and `abalones` downloads this morning (04:24–04:48 UTC) used the API, which is slower and requires pagination through 10k-record batches.
- S3 GeoParquet query is a single fast scan instead of repeated API calls; for `mammals_marine` (1.8M records), this saves significant time.
- AWS access verified working (360M total OBIS records accessible).
- API fallback ensures robustness if S3 is unavailable or returns nothing.

**Changes:**
- Added `config$obis_use_aws = TRUE` to `config/config.R` with documentation.
- Added `config$obis_max_age_days = 180` to match GBIF re-download age threshold.
- Updated `run_pipeline.R` to pass `use_aws = config$obis_use_aws` to `run_obis_downloads()`.
- Script already supported this via `download_one_taxon_obis(use_aws = TRUE)` default; just needed config to wire it in.

**Next run:** Expect OBIS downloads to show `source: "aws_geoparquet"` in manifests instead of `"api"`.

