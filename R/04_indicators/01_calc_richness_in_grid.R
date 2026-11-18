### =============================================================================
### R/04_indicators/01_calc_richness_in_grid.R  —  species richness per hex cell
### =============================================================================
###
### Summarise cleaned occurrences to per-hex-cell species richness and record
### counts, joined onto the region hex grid so the result carries geometry.
### Written as a GPKG per taxon-year (the sampling step reads it back).
###
### Caching: taxon-years whose clean parquets, hex grid and code are unchanged (and
### whose output still exists) are skipped; `force = TRUE` recomputes. Progress:
### `progress = "terse"` (default) one line per taxon-year, "verbose" adds detail,
### "quiet" prints nothing. Returns a status tibble (one row per taxon-year).
###
### Reworked from calc_gbif_spp_rich_in_grid.R: uses the canonical `h3_index`
### cell column, joins to the grid geometry (the original wrote a geometry-less
### tibble via st_write), and derives paths from sampling_config.
###
### Requires config + sampling_config + R/00_utils.R.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(purrr)
  library(sf)
  library(here)
  library(jsonlite)
})

### --- Cache -------------------------------------------------------------------
### A taxon-year is skipped when the md5 of its clean parquets, the hex grid and the
### function code all match the manifest and the output GPKG still exists.
### The manifest is written last, so an interrupted run is never treated as cached.
### `force = TRUE` recomputes everything.
richness_manifest_path <- function(taxon, year, cfg = config) {
  file.path(taxon_output_dir(taxon, cfg), "logs",
            glue::glue("richness_manifest_{taxon}_{year}.json"))
}

.read_richness_manifest <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

### One taxon-year: read every *_clean.parquet, summarise, join to grid, write.
### Returns a one-row status tibble (status: "written", "cached", "no_input").
calc_richness_in_grid <- function(taxon, year, cfg = config, samp_cfg = NULL,
                                  force = config$force,
                                  progress = c("quiet", "terse", "verbose")) {
  progress <- match.arg(progress)
  if (is.null(samp_cfg)) samp_cfg <- make_sampling_config(year, cfg = cfg)[[taxon]]
  occ_dir <- taxon_year_occ_dir(taxon, year, cfg)
  out_path <- samp_cfg$hex_richness_path
  status_row <- function(status, n_cells = NA_integer_) {
    tibble::tibble(taxon = taxon, year = year, status = status,
                   n_cells_with_spp = n_cells, file_path = out_path)
  }

  files <- list.files(occ_dir, pattern = "_clean\\.parquet$", full.names = TRUE)
  if (length(files) == 0) {
    warning(glue::glue("[{taxon} {year}] no _clean.parquet files in {occ_dir}"))
    return(status_row("no_input"))
  }

  fns <- c("calc_richness_in_grid")
  settings <- rlang::hash(list(
    inputs   = unname(tools::md5sum(files)),
    hex_grid = unname(tools::md5sum(hex_grid_path(cfg))),
    code     = lapply(fns, function(f) deparse(get(f, envir = globalenv())))
  ))
  mpath <- richness_manifest_path(taxon, year, cfg)
  old <- if (force) NULL else .read_richness_manifest(mpath)
  if (!is.null(old) && identical(old$settings, settings) && file.exists(out_path)) {
    .progress_msg(progress, "terse",
                  glue::glue("[richness] {taxon} {year}: CACHED (inputs unchanged; force = TRUE to rebuild)"))
    return(status_row("cached", old$n_cells_with_spp))
  }

  .progress_msg(progress, "terse",
                glue::glue("[richness] {taxon} {year} ({length(files)} files)"))
  occ <- purrr::map_dfr(files, ~arrow::read_parquet(
    .x, col_select = c("h3_index", "species", "tot_obs")))
  .progress_msg(progress, "verbose", "    read ", nrow(occ), " records")

  summary <- occ |>
    dplyr::filter(!is.na(species), !is.na(h3_index)) |>
    dplyr::group_by(h3_index) |>
    dplyr::summarise(n_species = dplyr::n_distinct(species),
                     n_records = sum(tot_obs, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(h3_index = as.character(h3_index))
  .progress_msg(progress, "verbose", "    ", nrow(summary), " cells with records")

  grid <- sf::st_read(hex_grid_path(cfg), quiet = TRUE)
  richness_sf <- dplyr::left_join(grid, summary, by = "h3_index") |>
    dplyr::mutate(n_species = tidyr::replace_na(n_species, 0L),
                  n_records = tidyr::replace_na(n_records, 0L))

  ensure_parent_dir(out_path)
  sf::st_write(richness_sf, out_path, delete_dsn = TRUE, quiet = TRUE)
  .progress_msg(progress, "verbose", "    wrote ", out_path)

  ### Manifest last.
  ensure_parent_dir(mpath)
  jsonlite::write_json(
    list(taxon = taxon, year = year, settings = settings,
         n_cells_with_spp = nrow(summary)),
    mpath, auto_unbox = TRUE, pretty = TRUE)
  status_row("written", nrow(summary))
}

run_richness_in_grid <- function(taxa = config$taxa, years = config$years, cfg = config,
                                 force = config$force,
                                 progress = c("terse", "verbose", "quiet")) {
  progress <- match.arg(progress)
  grid <- tidyr::expand_grid(taxon = taxa, year = years)
  purrr::map2_dfr(grid$taxon, grid$year,
                  ~calc_richness_in_grid(.x, .y, cfg, force = force, progress = progress))
}
