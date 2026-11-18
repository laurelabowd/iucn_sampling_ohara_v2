# h3o vs h3jsr benchmark for `range_from_polygons`

## Outcome

`h3o` (Rust H3 library) replaced `h3jsr` for the polygon-to-H3-cells step on 2026-11-18.
`R/03_process/04a_range_from_polygons_fn.R` now uses h3o; the previous h3jsr version is kept
at `R/03_process/_archive/04a_range_from_polygons_h3jsr_fn.R`.

Single benchmark pass, abalones, H3 resolution 4, centroid containment:

| | h3jsr | h3o |
|---|---|---|
| Total time (55 rows*) | 76.2 s | 5.0 s |
| Mean per species | 1.38 s | 0.09 s |
| Overall speedup | | 15.1x |
| Median per-species speedup | | 11.9x (IQR 9.9-14.8x) |
| Cell sets identical | | yes, all species |

\* The pre-scan found 43 species with >=100 cells, but a many-to-many join in the benchmark
script duplicates species with more than one crosswalk row (e.g. *Haliotis tuberculata*),
giving 55 rows. Timings are slightly inflated by the repeats; the ratio is unaffected in
practice. One run only (`N_RUNS <- 1`), so no variance estimate.

Raw results: `logs/benchmark_h3o_vs_h3jsr.csv`, `logs/benchmark_summary.txt`.
`logs/BENCHMARK_REPORT.md` and `logs/BENCHMARK_SUMMARY_TABLE.txt` are from before h3o was
installed and are stale (they report h3jsr only).

## API notes (h3o 0.3.0)

- h3o has **no** `polygon_to_cells()`. Use `h3o::sfc_to_cells(sf::st_geometry(x), resolution, containment)`
  then `h3o::flatten_h3()`.
- `containment = "centroid"` matches h3jsr (cells whose centre is inside the polygon). h3o's
  default, `"intersect"`, returns more cells and is not comparable.
- Input must be an `sfc_POLYGON`/`sfc_MULTIPOLYGON` in a degree CRS (EPSG:4326 here).
- Multi-feature species (e.g. breeding + foraging polygons) are converted together and
  flattened to one unique cell set, so each species is still one presence column.
- The strip-cutting fallback for antimeridian/invalid polygons was kept unchanged.

## Files

| File | Purpose |
|------|---------|
| `04a_range_from_polygons_fn.R` | Production range function (h3o) |
| `04_range_from_spatial_run.R` | Driver; sources 04a, unchanged function names |
| `05_benchmark_h3o_vs_h3jsr.R` | Benchmark harness (still needs both h3jsr and h3o) |
| `_archive/04a_range_from_polygons_h3jsr_fn.R` | Retired h3jsr implementation |
| `_archive/04c_range_from_polygons_h3o_fn.R` | Retired experimental h3o prototype (had no output-writing code) |

## Re-running the benchmark

```r
source("R/00_utils.R")
source("config/config.R")
source("config/taxa_config.R")
source("config/sampling_config.R")
source("R/03_process/05_benchmark_h3o_vs_h3jsr.R")
```

Settings at the top of the script: `TAXON`, `MIN_CELLS`, `N_RUNS` (use 3+ for variance),
`SAVE_CELLS`.

## Caching consequence

`04_range_from_spatial_run.R` hashes the deparsed range functions, so the code change
invalidates every polygon taxon's range manifest. The next `run_build_ranges()` rebuilds
them (h3o output should be cell-identical to the old output).

## Known issues

- Many-to-many join in the benchmark script (see above).
- h3o is a compiled package (needs Rust at install time); `renv.lock` should be updated
  with `renv::snapshot()` so the environment can be rebuilt.
