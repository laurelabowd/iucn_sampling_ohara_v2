### =============================================================================
### R/03_process/06_benchmark_hex_grid_h3o_vs_h3jsr.R  —  hex-grid build: h3o vs h3jsr
### =============================================================================
###
### Times the H3 part of build_hex_grid() (cell ids + cell polygons) with h3jsr and
### with h3o, and checks the two give identical cells and identical polygon vertices.
### The s2 conservation-area intersection is the same in both and is NOT timed.
### Writes nothing to the real grid layers.
###
### Input:  none (global grid built from res-0 cells; plus one regional bbox check)
### Output: logs/benchmark_hex_grid_h3o_vs_h3jsr.csv
###
### Usage (project root):  source("R/03_process/06_benchmark_hex_grid_h3o_vs_h3jsr.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(tibble)
  library(readr)
})

### 🔧 Settings ----
RES        <- 3     ### H3 level to benchmark
N_RUNS     <- 3     ### repeats per method (timings are wall-clock seconds)
### Regional-path correctness check only (not a data source): lon/lat bbox in degrees
TEST_BBOX  <- c(xmin = -130, xmax = -100, ymin = 20, ymax = 50)
OUT_CSV    <- here::here("logs", "benchmark_hex_grid_h3o_vs_h3jsr.csv")

### Vertex fingerprint per cell: rounded, sorted ring vertices (ignores start vertex)
.vertex_sets <- function(g) {
  lapply(sf::st_geometry(g), function(p) {
    m <- unclass(p)[[1]]
    m <- round(m[-nrow(m), , drop = FALSE], 6)
    m[order(m[, 1], m[, 2]), , drop = FALSE]
  })
}

### 🌍 1. Global grid ----
res0_h3o <- sprintf("8%03xfffffffffff", 2L * (0:121) + 1L)
stopifnot(identical(res0_h3o, h3jsr::get_res0()))   ### res-0 id formula matches h3jsr

time_global <- function(lib) {
  system.time({
    if (lib == "h3jsr") {
      ids <- h3jsr::get_children(h3jsr::get_res0(), res = RES, simple = TRUE) |>
        unlist(use.names = FALSE) |> unique()
      grid <- h3jsr::cell_to_polygon(ids, simple = FALSE) |>
        dplyr::rename(h3_index = h3_address)
    } else {
      ids <- res0_h3o |>
        h3o::h3_from_strings() |>
        h3o::get_children(resolution = RES) |>
        h3o::flatten_h3() |> as.character() |> unique()
      grid <- sf::st_sf(h3_index = ids, geometry = sf::st_as_sfc(h3o::h3_from_strings(ids)))
    }
  })[["elapsed"]]
}

timings <- tidyr::expand_grid(run = seq_len(N_RUNS), lib = c("h3jsr", "h3o")) |>
  dplyr::mutate(seconds = purrr::map_dbl(lib, time_global))

### One untimed build of each for the equality checks
ids_j  <- h3jsr::get_children(h3jsr::get_res0(), res = RES, simple = TRUE) |>
  unlist(use.names = FALSE) |> unique()
grid_j <- h3jsr::cell_to_polygon(ids_j, simple = FALSE) |> dplyr::rename(h3_index = h3_address)
ids_o  <- res0_h3o |> h3o::h3_from_strings() |> h3o::get_children(resolution = RES) |>
  h3o::flatten_h3() |> as.character() |> unique()
grid_o <- sf::st_sf(h3_index = ids_o, geometry = sf::st_as_sfc(h3o::h3_from_strings(ids_o)))

ids_same  <- setequal(ids_j, ids_o)
grid_o    <- grid_o[match(grid_j$h3_index, grid_o$h3_index), ]
vs_j <- .vertex_sets(grid_j); vs_o <- .vertex_sets(grid_o)
n_poly_same <- sum(mapply(function(x, y) identical(dim(x), dim(y)) &&
                            isTRUE(all.equal(x, y, check.attributes = FALSE)), vs_j, vs_o))

### 🗺️ 2. Regional path (polygon fill) correctness check ----
bb_poly <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(TEST_BBOX), crs = 4326)))
reg_j <- h3jsr::polygon_to_cells(bb_poly, res = RES, simple = TRUE) |> unlist(use.names = FALSE) |> unique()
reg_o <- h3o::sfc_to_cells(sf::st_geometry(bb_poly), resolution = RES, containment = "centroid") |>
  h3o::flatten_h3() |> as.character() |> unique()
reg_same <- setequal(reg_j, reg_o)

### 📋 3. Report ----
summary_tbl <- timings |>
  dplyr::summarise(n_runs = dplyr::n(), mean_sec = mean(seconds), min_sec = min(seconds),
                   max_sec = max(seconds), .by = lib)
speedup <- summary_tbl$mean_sec[summary_tbl$lib == "h3jsr"] / summary_tbl$mean_sec[summary_tbl$lib == "h3o"]

print(summary_tbl)
message(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "America/Los_Angeles"),
        " | level ", RES, ": ", length(ids_j), " cells (h3jsr) vs ", length(ids_o), " (h3o)")
message("Mean speedup: ", round(speedup, 1), "x")
message("Global cell ids identical: ", ids_same)
message("Cell polygons identical: ", n_poly_same, " / ", nrow(grid_j))
message("Regional bbox cell ids identical: ", reg_same, " (", length(reg_j), " vs ", length(reg_o), ")")

readr::write_csv(
  dplyr::mutate(timings, level = RES, n_cells = length(ids_j), ids_identical = ids_same,
                polygons_identical = n_poly_same == nrow(grid_j), regional_identical = reg_same),
  OUT_CSV)
