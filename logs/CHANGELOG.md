# CHANGELOG

## 2026-11-17 21:10 UTC — Consolidate OBIS download implementations

**Archived:** `R/02_acquire/02_download_obis_v2.R` and `02_download_obis_original_api_only.R` → `R/02_acquire/_archive/`

**Reason:** Consolidated to single canonical `02_download_obis.R` which supports both AWS S3 GeoParquet (primary) and API fallback. The archived versions were intermediate prototypes; the current version subsumes their functionality.

