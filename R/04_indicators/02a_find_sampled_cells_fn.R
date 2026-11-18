### =============================================================================
### R/04_indicators/02a_find_sampled_cells_fn.R  —  sampled-cell summaries (functions)
### =============================================================================
###
### For each species: how many of its range cells are sampled, how many of those
### fall inside conservation areas, and in/out-of-range observation counts.
###
### Reworked from find_sampled_cells_fn.R: the 57-line commented-out earlier copy
### is gone; everything keys off the canonical `h3_index` column; the threshold
### comes from sampling_params (one definition, one operator).
###
### Requires config + sampling_config + R/00_utils.R.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(purrr)
  library(readr)
  library(sf)
  library(here)
})

### One species file -> one-row summary. `cell_list` is the wide range-cell table
### (h3_index + one 0/1 column per species).
find_sampled_cells_one <- function(file_path, year, cons_ids, well_samp_cells, cell_list) {
  spp_df <- arrow::read_parquet(file_path,
                                col_select = c("h3_index", "species", "tot_obs"))
  spp_name <- unique(dplyr::pull(spp_df, species))

  na_row <- function(sp) tibble::tibble(
    species = sp, year = as.character(year),
    n_range_cells = NA_real_, tot_samp_cells = NA_real_, tot_samp_in_cons = NA_real_,
    n_obs_in_range = NA_real_, n_obs_out_of_range = NA_real_)

  if (length(spp_name) == 0) return(na_row(NA_character_))
  spp_name <- spp_name[[1]]
  if (!spp_name %in% names(cell_list)) return(na_row(spp_name))

  cell_list |>
    dplyr::select(h3_index, dplyr::all_of(spp_name)) |>
    dplyr::left_join(spp_df, by = "h3_index") |>
    dplyr::mutate(
      in_range  = .data[[spp_name]] == 1 & !is.na(.data[[spp_name]]),
      has_obs   = !is.na(tot_obs),
      cons_area = h3_index %in% cons_ids,
      well_samp = h3_index %in% well_samp_cells,
      samp_cell = well_samp | has_obs
    ) |>
    dplyr::summarize(
      species = spp_name, year = as.character(year),
      n_range_cells = sum(in_range, na.rm = TRUE),
      tot_samp_cells = sum(samp_cell, na.rm = TRUE),
      tot_samp_in_cons = sum(samp_cell & cons_area, na.rm = TRUE),
      n_obs_in_range = sum(has_obs & in_range, na.rm = TRUE),
      n_obs_out_of_range = sum(has_obs & !in_range, na.rm = TRUE)
    )
}

### Full workflow for one taxon-year.
### Cached on (clean parquets, range cells, conservation hexes, richness gpkg,
### threshold, code); a cached run re-reads the CSV. `force = TRUE` recomputes.
### `progress`: "quiet" | "terse" | "verbose".
run_find_sampled_cells <- function(taxon, year, cfg = config,
                                   params = sampling_params, write_output = TRUE,
                                   force = config$force,
                                   progress = c("quiet", "terse", "verbose")) {
  progress <- match.arg(progress)
  samp_cfg <- make_sampling_config(year, cfg = cfg)[[taxon]]
  out <- file.path(samp_cfg$output_dir, glue::glue("sampled_cells_{taxon}_{year}.csv"))

  occ_dir <- taxon_year_occ_dir(taxon, year, cfg)
  files <- list.files(occ_dir, pattern = "_clean\\.parquet$", full.names = TRUE)

  mpath <- file.path(taxon_output_dir(taxon, cfg), "logs",
                     glue::glue("sampled_cells_manifest_{taxon}_{year}.json"))
  fns <- c("find_sampled_cells_one", "run_find_sampled_cells")
  settings <- rlang::hash(list(
    clean  = unname(tools::md5sum(files)),
    others = unname(tools::md5sum(c(samp_cfg$range_cells_path, cons_areas_hex_path(cfg),
                                    samp_cfg$hex_richness_path))),
    well_sampled_threshold = params$well_sampled_threshold,
    code   = lapply(fns, function(f) deparse(get(f, envir = globalenv())))
  ))
  if (isTRUE(write_output) && .stage_cached(mpath, settings, out, force)) {
    .progress_msg(progress, "terse",
                  glue::glue("[sampled cells] {taxon} {year}: CACHED (inputs unchanged; force = TRUE to rebuild)"))
    return(readr::read_csv(out, show_col_types = FALSE))
  }

  .progress_msg(progress, "terse",
                glue::glue("[sampled cells] {taxon} {year} ({length(files)} species files)"))
  cell_list <- arrow::read_parquet(samp_cfg$range_cells_path)
  cons_ids  <- sf::st_read(cons_areas_hex_path(cfg), quiet = TRUE)$h3_index

  richness <- sf::st_read(samp_cfg$hex_richness_path, quiet = TRUE) |> sf::st_drop_geometry()
  well_samp_cells <- richness |>
    dplyr::filter(n_records >= params$well_sampled_threshold) |>
    dplyr::pull(h3_index)
  .progress_msg(progress, "verbose", "    ", length(well_samp_cells), " well-sampled cells")

  results <- purrr::map_dfr(files, find_sampled_cells_one,
                            year = year, cons_ids = cons_ids,
                            well_samp_cells = well_samp_cells, cell_list = cell_list)
  .progress_msg(progress, "verbose", "    summarised ", nrow(results), " species")

  if (isTRUE(write_output)) {
    ensure_parent_dir(out)
    readr::write_csv(results, out)
    .write_stage_manifest(mpath, settings, taxon = taxon, year = year)
    .progress_msg(progress, "verbose", "    wrote ", out)
  }
  results
}
