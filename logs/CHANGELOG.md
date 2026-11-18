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
