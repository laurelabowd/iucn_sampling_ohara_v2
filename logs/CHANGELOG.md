# CHANGELOG

## 2026-11-17 21:10 UTC — Consolidate OBIS download implementations

**Archived:** `R/02_acquire/02_download_obis_v2.R` and `02_download_obis_original_api_only.R` → `R/02_acquire/_archive/`

**Reason:** Consolidated to single canonical `02_download_obis.R` which supports both AWS S3 GeoParquet (primary) and API fallback. The archived versions were intermediate prototypes; the current version subsumes their functionality.


## 2026-11-18 PST — Switch range_from_polygons to h3o

**Archived:** `R/03_process/04a_range_from_polygons_fn.R` (h3jsr version) → `R/03_process/_archive/04a_range_from_polygons_h3jsr_fn.R`
**Archived:** `R/03_process/04c_range_from_polygons_h3o_fn.R` (experimental prototype) → `R/03_process/_archive/`
**Edited:** `R/03_process/04a_range_from_polygons_fn.R` now uses `h3o::sfc_to_cells()`; function names unchanged, so `04_range_from_spatial_run.R` needed no edit.
**Edited:** `R/03_process/05_benchmark_h3o_vs_h3jsr.R` (h3o call fixed to `sfc_to_cells`), `R/03_process/README_h3o_benchmark.md` (rewritten with results).
**Note:** range manifests hash the range code, so all polygon taxa will rebuild on the next run.

## 2026-11-18 PST — Switch build_hex_grid to h3o; add hex-grid benchmark

**Edited:** `R/03_process/03_build_hex_grid.R` — h3jsr calls (`get_children`/`get_res0`, `polygon_to_cells`, `cell_to_polygon`) replaced with h3o (`get_children`, `sfc_to_cells`, `sf::st_as_sfc()` on H3 vectors).
**Added:** `R/03_process/06_benchmark_hex_grid_h3o_vs_h3jsr.R` → `logs/benchmark_hex_grid_h3o_vs_h3jsr.csv`.
**Note:** the hex-grid manifest hashes inputs, not code, so existing grid layers stay cached; use `build_hex_grid(force = TRUE)` to regenerate with h3o (not yet done).

## 2026-11-18 09:37 PST — Add caching and progress option to richness-in-grid

**Edited:** `R/04_indicators/01_calc_richness_in_grid.R` — `calc_richness_in_grid()` / `run_richness_in_grid()` gain `force` and `progress` ("terse"/"verbose"/"quiet") args. Per taxon-year manifest `<taxon output dir>/logs/richness_manifest_<taxon>_<year>.json` (hash of clean-parquet md5s, hex-grid md5, function code). Both functions now return a status tibble instead of the sf / `TRUE`.
**Edited:** `R/run_pipeline.R` — call is now `run_richness_in_grid(force = FALSE)`; stale TODO note removed.
**Edited:** `R/README.md` — richness row notes caching/progress.
**Note:** no existing manifests, so the first run after this change recomputes every taxon-year.

## 2026-11-18 10:15 PST — Add caching and progress option to sampled-cells, optimal-locs and greedy set-cover

**Edited:** `R/00_utils.R` — new shared helpers `.stage_cached()` and `.write_stage_manifest()` (manifest = settings hash; reusable only if hash matches and outputs exist; manifest written last).
**Edited:** `R/04_indicators/02a_find_sampled_cells_fn.R`, `02_find_sampled_cells_run.R`, `03a_find_optimal_sampling_locs_fn.R`, `03_find_optimal_sampling_locs_run.R`, `04a_greedy_setcover_fn.R`, `04_greedy_setcover_run.R` — `force` and `progress` ("terse"/"verbose"/"quiet") args. Manifests: `<taxon output dir>/logs/{sampled_cells,optimal_locs}_manifest_<taxon>_<year>.json`; `<greedy output dir>/greedy_manifest.json`. Hashes cover upstream input md5s, relevant `sampling_params`, and function code.
**Edited:** `R/run_pipeline.R` — stage 3 calls pass `force = FALSE`.
**Note:** `find_optimal_sampling_locs()` and `greedy_setcover()` now return a `status` element ("written"/"cached"); cached runs return only `status` plus paths, not the in-memory tables. No existing manifests, so the first run recomputes everything.

## 2026-11-18 10:40 PST — `force` defaults to `config$force` in all stage functions

**Edited:** every `force = FALSE` default in `R/01_taxonomy`, `R/02_acquire`, `R/03_process`, `R/04_indicators` (not the benchmark helpers or internal `.stage_cached()`) is now `force = config$force`. Explicit `force = ...` arguments in `R/run_pipeline.R` still override it.
