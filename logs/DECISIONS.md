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


## 2026-11-18 PST — Use h3o instead of h3jsr for polygon-to-H3 conversion

**Decision:** Production range step uses `h3o::sfc_to_cells(containment = "centroid")`.

**Evidence:** Abalones, res 4, one pass: h3jsr 76.2 s vs h3o 5.0 s (15.1x overall; median 11.9x per species, IQR 9.9-14.8x). Cell sets identical for every species benchmarked. Only one run, no variance estimate; only abalones tested.

**Rationale:** Same cells, an order of magnitude faster. `"centroid"` containment chosen because h3jsr keeps cells whose centre lies in the polygon; h3o's default `"intersect"` would change range sizes.

**Rejected / ruled out:** The prototype assumed `h3o::polygon_to_cells()`; it does not exist in h3o 0.3.0. Batching all species in one call (grouping cells by species afterwards) was discussed but not implemented; multi-feature species are still combined per species before flattening.

**Open:** Not yet verified on other polygon taxa, especially antimeridian-crossing ranges that triggered the strip-cutting fallback with h3jsr. Benchmark script has a many-to-many join inflating 43 species to 55 rows. `renv.lock` not yet re-snapshotted for h3o.
