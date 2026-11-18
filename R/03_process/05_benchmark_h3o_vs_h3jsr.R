### =============================================================================
### R/03_process/05_benchmark_h3o_vs_h3jsr.R  —  compare polygon-to-cells speed
### =============================================================================
###
### Benchmark h3o::polygon_to_cells against h3jsr::polygon_to_cells for the
### abalone taxon. Compares wall-clock time, cell counts, and cell identity
### across a sample of species.
###
### Input: IUCN polygon files (crosswalk, range shapefiles, config, taxa_config)
### Output: logs/benchmark_h3o_vs_h3jsr.csv, logs/benchmark_summary.txt
###
### Usage:
###   source("R/00_utils.R")
###   source("config/config.R")
###   source("config/taxa_config.R")
###   source("config/sampling_config.R")
###   source("R/03_process/05_benchmark_h3o_vs_h3jsr.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(glue)
  library(sf)
  library(h3jsr)
  library(tibble)
  library(here)
})

### Settings
TAXON        <- "abalones"
MIN_CELLS    <- 100          ### only benchmark species with >100 cells (meaningful runtime)
N_RUNS       <- 1            ### number of times to run each species (1 = single pass)
SAVE_CELLS   <- TRUE         ### save cell sets for comparison

### Ensure logs directory exists
ensure_dir("logs")

###  =============================================================================
### Stage 1: Load crosswalk and identify polygon-linked species
###  =============================================================================

.progress_msg("terse", "terse", glue("[{TAXON}] Loading crosswalk and polygon metadata..."))

cw <- readr::read_csv(config$paths$crosswalk, show_col_types = FALSE)
if (!"range_polygon_layer" %in% names(cw)) cw$range_polygon_layer <- NA_character_

spp_map <- cw |>
  dplyr::filter(taxon == TAXON, !is.na(range_polygon_name), !is.na(name_gbif)) |>
  dplyr::transmute(range_polygon_name, species = name_gbif,
                   polygon_path = range_polygon_path,
                   polygon_layer = as.character(range_polygon_layer)) |>
  tidyr::separate_longer_delim(c(range_polygon_name, polygon_path, polygon_layer),
                               delim = "|")

if (nrow(spp_map) == 0) {
  stop(glue("[{TAXON}] no polygon-linked species in crosswalk."), call. = FALSE)
}

tcfg <- taxa_config[[TAXON]]
n_sources <- dplyr::distinct(spp_map, polygon_path, polygon_layer) |> nrow()
n_files   <- dplyr::n_distinct(spp_map$polygon_path)

.progress_msg("terse", "terse",
              glue("[{TAXON}] found {dplyr::n_distinct(spp_map$species)} species ",
                   "from {n_files} file{if (n_files != 1) 's' else ''}",
                   "{if (n_sources > n_files) glue(' / {n_sources} layers') else ''}"))

###  =============================================================================
### Stage 2: Pre-scan with h3jsr to identify species with enough cells
###  =============================================================================

.progress_msg("terse", "terse", glue("[{TAXON}] Pre-scanning with h3jsr to filter by cell count..."))

### Read and process each layer once
species_cell_counts <- spp_map |>
  dplyr::group_split(polygon_path, polygon_layer, .keep = TRUE) |>
  purrr::map_dfr(function(g) {
    path <- g$polygon_path[1]
    layer <- g$polygon_layer[1]
    lyr_lab <- if (is.na(layer)) "default" else layer
    
    tryCatch({
      range_sf <- if (is.na(layer)) sf::st_read(path, quiet = TRUE)
                  else sf::st_read(path, layer = layer, quiet = TRUE)
      range_sf <- sf::st_transform(range_sf, config$crs$working)
      
      purrr::pmap_dfr(dplyr::select(g, range_polygon_name, species),
                     function(range_polygon_name, species) {
        sp_sf <- range_sf[as.character(range_sf[[tcfg$polygon_name_field]]) %in% range_polygon_name, ]
        if (nrow(sp_sf) == 0) return(NULL)
        
        cells <- tryCatch(
          h3jsr::polygon_to_cells(sp_sf, res = config$hex_resolution) |> unlist(),
          error = function(e) character(0)
        )
        
        tibble::tibble(
          species = species,
          n_cells = length(unique(cells)),
          n_features = nrow(sp_sf),
          file = basename(path),
          layer = lyr_lab
        )
      })
    }, error = function(e) {
      warning(glue("failed to scan {basename(path)} [{lyr_lab}]: {conditionMessage(e)}"))
      NULL
    })
  }) |>
  dplyr::filter(n_cells >= MIN_CELLS) |>
  dplyr::arrange(dplyr::desc(n_cells))

.progress_msg("terse", "terse",
              glue("[{TAXON}] {nrow(species_cell_counts)} species have ≥{MIN_CELLS} cells"))
print(species_cell_counts, n = 15)

if (nrow(species_cell_counts) == 0) {
  stop(glue("[{TAXON}] no species with ≥{MIN_CELLS} cells; nothing to benchmark."),
       call. = FALSE)
}

###  =============================================================================
### Stage 3: Benchmark h3jsr and h3o on each qualifying species
###  =============================================================================

### Helper: convert one species' polygons using a specified library
.benchmark_one_species <- function(range_sf, range_polygon_name, species, lib = c("h3jsr", "h3o")) {
  lib <- match.arg(lib)
  
  sp_sf <- range_sf[as.character(range_sf[[tcfg$polygon_name_field]]) %in% range_polygon_name, ]
  if (nrow(sp_sf) == 0) return(list(cells = character(0), time_sec = NA_real_))
  
  t0 <- Sys.time()
  cells <- tryCatch({
    if (lib == "h3jsr") {
      suppressMessages(h3jsr::polygon_to_cells(sp_sf, res = config$hex_resolution)) |> unlist()
    } else {
      ### h3o has no polygon_to_cells(); the equivalent is sfc_to_cells(). h3jsr
      ### uses centroid containment, so match it for an apples-to-apples result.
      if (!requireNamespace("h3o", quietly = TRUE)) {
        return(list(cells = character(0), time_sec = NA_real_))
      }
      h3o::sfc_to_cells(sf::st_geometry(sp_sf), resolution = config$hex_resolution,
                        containment = "centroid") |>
        h3o::flatten_h3() |>
        as.character()
    }
  }, error = function(e) {
    warning(glue("[{species}] {lib} failed: {conditionMessage(e)}"))
    character(0)
  })
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(cells = as.character(unique(cells)), time_sec = elapsed)
}

### Benchmark on each layer, once per library
benchmark_results <- list()

for (run in seq_len(N_RUNS)) {
  .progress_msg("terse", "terse", glue("[run {run}/{N_RUNS}] Reading layers and benchmarking..."))
  
  run_results <- spp_map |>
    dplyr::group_split(polygon_path, polygon_layer, .keep = TRUE) |>
    purrr::map_dfr(function(g) {
      path <- g$polygon_path[1]
      layer <- g$polygon_layer[1]
      
      ### Only process species we identified in pre-scan
      g <- g |> dplyr::inner_join(
        species_cell_counts |> dplyr::select(species),
        by = "species"
      )
      if (nrow(g) == 0) return(NULL)
      
      tryCatch({
        range_sf <- if (is.na(layer)) sf::st_read(path, quiet = TRUE)
                    else sf::st_read(path, layer = layer, quiet = TRUE)
        range_sf <- sf::st_transform(range_sf, config$crs$working)
        
        purrr::pmap_dfr(dplyr::select(g, range_polygon_name, species),
                       function(range_polygon_name, species) {
          result_h3jsr <- .benchmark_one_species(range_sf, range_polygon_name, species, lib = "h3jsr")
          result_h3o <- .benchmark_one_species(range_sf, range_polygon_name, species, lib = "h3o")
          
          tibble::tibble(
            run = run,
            species = species,
            n_cells_h3jsr = length(result_h3jsr$cells),
            n_cells_h3o = length(result_h3o$cells),
            time_h3jsr_sec = result_h3jsr$time_sec,
            time_h3o_sec = result_h3o$time_sec,
            speedup = if (!is.na(result_h3jsr$time_sec) && !is.na(result_h3o$time_sec)) {
              result_h3jsr$time_sec / result_h3o$time_sec
            } else NA_real_
          ) |>
            dplyr::mutate(
              cells_match = if (length(result_h3jsr$cells) > 0 && length(result_h3o$cells) > 0) {
                length(intersect(result_h3jsr$cells, result_h3o$cells)) == length(result_h3jsr$cells) &&
                length(intersect(result_h3jsr$cells, result_h3o$cells)) == length(result_h3o$cells)
              } else NA
            )
        })
      }, error = function(e) {
        warning(glue("failed to process file {basename(path)}: {conditionMessage(e)}"))
        NULL
      })
    })
  
  benchmark_results[[run]] <- run_results
}

benchmark_df <- purrr::list_rbind(benchmark_results)

###  =============================================================================
### Stage 4: Summarize and save results
###  =============================================================================

.progress_msg("terse", "terse", glue("[{TAXON}] Writing results..."))

### Detailed results to CSV
readr::write_csv(
  benchmark_df |>
    dplyr::mutate(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      h3o_available = !all(is.na(time_h3o_sec))
    ) |>
    dplyr::select(timestamp, species, n_cells_h3jsr, n_cells_h3o, 
                  time_h3jsr_sec, time_h3o_sec, speedup, cells_match),
  "logs/benchmark_h3o_vs_h3jsr.csv"
)

### Summary statistics
summary_stats <- benchmark_df |>
  dplyr::summarise(
    n_species = dplyr::n(),
    h3jsr_total_time = sum(time_h3jsr_sec, na.rm = TRUE),
    h3o_total_time = sum(time_h3o_sec, na.rm = TRUE),
    h3jsr_mean_time = mean(time_h3jsr_sec, na.rm = TRUE),
    h3o_mean_time = mean(time_h3o_sec, na.rm = TRUE),
    median_speedup = median(speedup, na.rm = TRUE),
    q25_speedup = quantile(speedup, 0.25, na.rm = TRUE),
    q75_speedup = quantile(speedup, 0.75, na.rm = TRUE),
    cells_match_all = all(cells_match == TRUE, na.rm = TRUE),
    .by = run
  )

### Write summary
summary_file <- "logs/benchmark_summary.txt"
writeLines(
  c(
    "================================================================================",
    glue("H3O vs H3JSR Benchmark: {TAXON}"),
    glue("Run at: {format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z', tz = 'America/Los_Angeles')}"),
    "================================================================================",
    "",
    glue("Settings:"),
    glue("  Taxon: {TAXON}"),
    glue("  Min cells per species: {MIN_CELLS}"),
    glue("  H3 resolution: {config$hex_resolution}"),
    glue("  Number of runs: {N_RUNS}"),
    "",
    glue("Results ({nrow(summary_stats)} run{if (nrow(summary_stats) > 1) 's' else ''}):"),
    "",
    purrr::pmap_chr(summary_stats, function(run, n_species, h3jsr_total_time, 
                                           h3o_total_time, h3jsr_mean_time,
                                           h3o_mean_time, median_speedup,
                                           q25_speedup, q75_speedup, cells_match_all) {
      h3o_available <- !is.na(h3o_total_time)
      if (!h3o_available) {
        return(glue("  Run {run}: {n_species} species | h3jsr: {format(round(h3jsr_total_time, 1), nsmall=1)}s | h3o: NOT INSTALLED"))
      }
      overall_speedup <- h3jsr_total_time / h3o_total_time
      paste(c(
        glue("  Run {run}:"),
        glue("    Species benchmarked: {n_species}"),
        glue("    h3jsr total: {format(round(h3jsr_total_time, 1), nsmall=1)}s (mean: {format(round(h3jsr_mean_time, 2), nsmall=2)}s/spp)"),
        glue("    h3o total:   {format(round(h3o_total_time, 1), nsmall=1)}s (mean: {format(round(h3o_mean_time, 2), nsmall=2)}s/spp)"),
        glue("    Overall speedup: {format(round(overall_speedup, 2), nsmall=2)}x"),
        glue("    Median speedup per species: {format(round(median_speedup, 2), nsmall=2)}x (IQR: {format(round(q25_speedup, 2), nsmall=2)} – {format(round(q75_speedup, 2), nsmall=2)}x)"),
        glue("    Cell sets match: {if (cells_match_all) 'YES ✓' else 'NO ✗'}")
      ), collapse = "\n")
    }),
    "",
    "Details by species: logs/benchmark_h3o_vs_h3jsr.csv",
    "================================================================================",
    ""
  ),
  summary_file
)

cat(readr::read_file(summary_file))

### Write to project memory
.progress_msg("terse", "terse", glue("Benchmark complete. Summary in {summary_file}"))
