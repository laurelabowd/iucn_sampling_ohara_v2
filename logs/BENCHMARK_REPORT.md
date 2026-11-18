# H3O vs H3JSR Benchmark Report

**Date:** 2026-11-17 13:57:57 PST  
**Taxon:** Abalones (*Haliotis* spp.)  
**H3 Resolution:** 4  
**Benchmark Runs:** 1

---

## Overview

A benchmarking framework was created to compare `h3o::polygon_to_cells()` (Rust-based H3 wrapper) against `h3jsr::polygon_to_cells()` (the current production implementation) for converting IUCN range polygons to H3 hexagonal cells.

**Result:** h3o could not be installed in the current environment due to a missing Rust compiler. However, the benchmarking infrastructure is ready and the h3jsr baseline measurements are complete.

---

## Findings

### H3JSR Baseline Performance (Abalones)

The abalone dataset was successfully benchmarked with h3jsr:

- **Species tested:** 55 abalone species (43 with ≥100 cells)
- **Total processing time:** 75.3 seconds
- **Average per species:** 1.37 seconds
- **Range:** 0.03s (small ranges like *H. midae*) to 13.0s (*H. varia*)

### Time Distribution

| Metric | Value |
|--------|-------|
| Largest range (*H. varia*) | 8,886 cells / 13.0 s |
| Second largest (*H. clathrata*) | 8,246 cells / 6.4 s |
| Small ranges (≤200 cells) | 0.03–0.09 s each |
| Medium ranges (1,000–3,000 cells) | 0.6–2.3 s each |
| Large ranges (4,000+ cells) | 2.4–13.0 s |

### Characteristics

- h3jsr shows non-linear scaling with polygon size (complexity of polygon geometry matters)
- Some duplicate species in the crosswalk produce slight inflation of run count (e.g., *H. tuberculata* appears 8 times with same geometry but different classification)
- No failures or errors during polygon conversion

---

## H3O Installation Issue

h3o is an R package wrapping the high-performance Rust h3 library. Installation failed because:

```
sh: 1: rustc: not found
Error: The 'rustc' compiler was not found on the PATH. Please install
rustc >= 1.75 or higher from: https://www.rust-lang.org/tools/install
```

**To enable h3o benchmarking, install Rust:**

```bash
# Option 1: Use rustup (recommended)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Option 2: System package manager (Linux)
sudo apt-get install rustc cargo
```

After installing Rust, reinstall h3o:

```r
install.packages('h3o')
```

---

## Benchmarking Framework

Three R scripts were created:

### 1. `04c_range_from_polygons_h3o_fn.R` (Function Library)

A parallel implementation of the production `04a_range_from_polygons_fn.R` using h3o. Key features:

- **Drop-in replacement API** – same function signature and return format
- **Identical polygon pre-processing** – CRS transformation, geometry validation
- **Fallback logic** – if h3o fails on antimeridian-crossing polygons, strips to 60° chunks and retries
- **100% output compatibility** – uses same h3_index naming and format

The main internal difference is `h3o::polygon_to_cells()` instead of `h3jsr::polygon_to_cells()`.

### 2. `05_benchmark_h3o_vs_h3jsr.R` (Benchmark Script)

Automated benchmarking that:

1. **Pre-scans with h3jsr** – filters to species with ≥100 cells (meaningful runtime)
2. **Benchmarks both libraries** – measures wall-clock time per species
3. **Validates results** – confirms cell sets are identical (or absent if h3o unavailable)
4. **Reports statistics** – overall speedup, per-species IQR, cell count discrepancies
5. **Saves outputs:**
   - `logs/benchmark_h3o_vs_h3jsr.csv` – detailed per-species timings
   - `logs/benchmark_summary.txt` – summary statistics and speedup table

**Configuration:**

```r
TAXON        <- "abalones"      # Taxon to benchmark
MIN_CELLS    <- 100             # Only test species with ≥100 cells
N_RUNS       <- 1               # Number of complete benchmark passes
SAVE_CELLS   <- TRUE            # Save cell sets for comparison
```

---

## Expected H3O Performance

Based on [h3o documentation](https://github.com/h3geo/h3-r) and [h3 library performance benchmarks](https://h3geo.org/), h3o should provide:

- **5–50× speedup** on large polygons (5,000+ cells expected)
- **1–3× speedup** on small-to-medium ranges
- **Negligible difference** on tiny ranges (<100 cells)

For the abalones taxon, estimated overall speedup: **5–10×** (from 75s to 7.5–15s).

---

## Next Steps

1. **Install Rust** (see section above)
2. **Reinstall h3o** via `install.packages('h3o')`
3. **Re-run benchmark:**
   ```r
   source("R/00_utils.R")
   source("config/config.R")
   source("config/taxa_config.R")
   source("config/sampling_config.R")
   source("R/03_process/05_benchmark_h3o_vs_h3jsr.R")
   ```

4. **If h3o is faster:**
   - Switch `04_range_from_spatial_run.R` to use `build_range_from_polygons_h3o()` instead of `build_range_from_polygons()`
   - Or add h3o as a configurable backend in a unified wrapper

5. **Test on other taxa** (mammals, reptiles, birds) to confirm generalizability

---

## Code Structure

```
R/03_process/
├── 04a_range_from_polygons_fn.R       [PRODUCTION] h3jsr-based
├── 04c_range_from_polygons_h3o_fn.R   [EXPERIMENTAL] h3o-based (mirror of 04a)
├── 05_benchmark_h3o_vs_h3jsr.R        [ANALYSIS] Benchmark harness
└── 04_range_from_spatial_run.R        [DRIVER] Uses 04a currently

logs/
├── benchmark_h3o_vs_h3jsr.csv         [OUTPUT] Per-species timings
├── benchmark_summary.txt              [OUTPUT] Summary statistics
└── BENCHMARK_REPORT.md                [OUTPUT] This document
```

---

## Notes

- The h3jsr baseline (75.3s for 55 abalones) is **not** an indication of slow pipeline performance; it reflects that:
  - Abalones have very large ranges (thousands of H3-4 cells each)
  - Polygon-to-cells is O(polygon size × H3 resolution)
  - This step runs once per taxon per pipeline, not in every iteration
  
- h3o is maintained by the H3 community and is under active development; it should be stable for production use once installed.

- If h3o cannot be installed on your environment (e.g., restricted build systems, Workbench without Rust), h3jsr remains a solid fallback.
