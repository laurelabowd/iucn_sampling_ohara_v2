### =============================================================================
### R/03_process/02_clean_occurrence.R  —  clean + spatially thin occurrences (GBIF and/or OBIS)
### =============================================================================
###
### For each per-species-year parquet from 01_parse_occurrence.R: run CoordinateCleaner,
### apply quality filters (thresholds from config$cleaning), diagnostic checks,
### and H3 thinning (one record per hex cell at config$hex_resolution). Writes
### `<species>_clean.parquet` next to the input and a per-taxon-year run log.
###
### Reworked from the original clean_gbif_occs.R: thresholds and CRS/resolution
### come from config; the ~13 hand-toggled input/output pairs are replaced by a
### loop over taxa x years; every species' outcome is recorded in a run log
### instead of a bare cat()/NULL.
###
### Performance: the country-centroid (cc_cen) and capital (cc_cap) tests are done by
### `cc_near_ref()` against geodesic buffers built once per R session; the
### CoordinateCleaner versions rebuild them for every file (~7 s each). Results are
### identical. Thinning works on plain columns (no sf round trip).
###
### Caching: per taxon-year, a manifest (clean_occurrence_manifest_<taxon>_<year>.json
### in the taxon's logs folder) records each input file's md5 plus a hash of
### config$cleaning, hex resolution, CRS and the cleaning code. Files whose input and
### settings are unchanged (and whose output still exists) are skipped. `force = TRUE`
### recleans everything.
###
### Progress messages: `progress = "terse"` (default) prints one line per taxon-year,
### `"verbose"` adds per-file and per-step lines, `"quiet"` prints nothing.
###
### Requires config + R/00_utils.R.  Usage:  res <- run_clean_occurrence()
###                                          res <- run_clean_occurrence(force = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(glue)
  library(arrow)
  library(sf)
  library(CoordinateCleaner)
  library(h3jsr)
  library(here)
})

### GBIF SIMPLE_PARQUET uses lower-case column names; normalise the ones we use to
### the canonical camelCase, renaming only those actually present (the original
### rename() errored if any were missing).
.normalise_gbif_cols <- function(df) {
  map <- c(
    gbifID = "gbifid", datasetKey = "datasetkey", occurrenceID = "occurrenceid",
    scientificName = "scientificname", occurrenceStatus = "occurrencestatus",
    individualCount = "individualcount", publishingOrgKey = "publishingorgkey",
    decimalLatitude = "decimallatitude", decimalLongitude = "decimallongitude",
    coordinateUncertaintyInMeters = "coordinateuncertaintyinmeters",
    coordinatePrecision = "coordinateprecision", eventDate = "eventdate",
    basisOfRecord = "basisofrecord"
  )
  for (canonical in names(map)) {
    lower <- map[[canonical]]
    if (!canonical %in% names(df) && lower %in% names(df)) {
      names(df)[names(df) == lower] <- canonical
    }
  }
  df
}

### --- Reference-buffer tests (cc_cen / cc_cap replacement) ---------------------
### Same test as CoordinateCleaner::cc_cen()/cc_cap() (geod = TRUE): flag points
### within `buffer` metres (geodesic) of a country centroid / capital. CoordinateCleaner
### rebuilds the buffer polygons on every call; here they are built once per buffer.
.wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"
.ref_buffer_cache <- new.env(parent = emptyenv())

.ref_buffers <- function(kind = c("centroid", "capital"), buffer) {
  kind <- match.arg(kind)
  key <- paste(kind, buffer)
  if (!is.null(.ref_buffer_cache[[key]])) return(.ref_buffer_cache[[key]])

  ref <- CoordinateCleaner::countryref
  cols <- if (kind == "centroid") c("centroid.lon", "centroid.lat") else c("capital.lon", "capital.lat")
  if (kind == "capital") ref <- ref[!is.na(ref$capital), ]
  ref <- ref[stats::complete.cases(ref[, cols]), cols]

  dg <- seq(from = 0, to = 360, by = 5)
  n <- nrow(ref)
  xy <- geosphere::destPoint(p = as.matrix(ref), b = rep(dg, each = n), d = buffer)
  id <- rep(seq_len(n), times = length(dg))
  ring <- order(id)
  polys <- terra::vect(cbind(id[ring], 1, xy[ring, 1], xy[ring, 2]),
                       type = "polygons", crs = .wgs84)
  .ref_buffer_cache[[key]] <- polys
  polys
}

### Drop records within `buffer` metres of a reference point (kind = "centroid" or
### "capital"). Drop-in for cc_cen(buffer = ) / cc_cap(buffer = ).
cc_near_ref <- function(d, kind, buffer, lon = "decimalLongitude", lat = "decimalLatitude") {
  if (nrow(d) == 0) return(d)
  pts <- terra::vect(as.data.frame(d[, c(lon, lat)]), geom = c(lon, lat), crs = .wgs84)
  hit <- terra::extract(.ref_buffers(kind, buffer), pts)
  d[is.na(hit[!duplicated(hit[, 1]), 2]), ]
}

### Print `...` if `progress` is at least `level` ("terse" < "verbose").
.progress_msg <- function(progress, level, ...) {
  rank <- c(quiet = 0, terse = 1, verbose = 2)
  if (rank[[progress]] >= rank[[level]]) message(...)
  invisible()
}

### Clean + thin one species file. Returns a one-row status tibble (never stops).
clean_one_file <- function(input_path, output_dir, cfg = config,
                           progress = c("quiet", "terse", "verbose")) {
  progress <- match.arg(progress)
  clean_name <- basename(input_path) |>
    str_remove("_occs\\.parquet$") |>
    janitor::make_clean_names()

  res <- tryCatch({
    occ <- arrow::read_parquet(input_path) |> .normalise_gbif_cols()
    ### Files from before OBIS merging have no record_id; fall back to gbifID.
    if (!"record_id" %in% names(occ)) occ$record_id <- as.character(occ$gbifID)
    cl  <- cfg$cleaning
    .progress_msg(progress, "verbose", "    read ", nrow(occ), " records")

    d <- occ |>
      dplyr::filter(!basisOfRecord %in% cl$exclude_basis_of_record) |>
      dplyr::filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) |>
      cc_dupl(verbose = FALSE) |> cc_zero(verbose = FALSE) |>
      cc_equ(verbose = FALSE) |> cc_val(verbose = FALSE) |>
      cc_near_ref("capital", buffer = 2000) |>
      cc_near_ref("centroid", buffer = 2000) |>
      cc_gbif(buffer = 2000, verbose = FALSE) |> cc_inst(buffer = 500, verbose = FALSE) |>
      dplyr::filter(coordinateUncertaintyInMeters < cl$max_coord_uncertainty_m |
                      is.na(coordinateUncertaintyInMeters)) |>
      dplyr::mutate(coordinatePrecision = as.numeric(coordinatePrecision)) |>
      dplyr::filter(coordinatePrecision > cl$min_coord_precision |
                      is.na(coordinatePrecision)) |>
      dplyr::mutate(individualCount = as.numeric(individualCount)) |>
      dplyr::filter(individualCount > 0 | is.na(individualCount)) |>
      dplyr::filter(individualCount < cl$max_individual_count | is.na(individualCount))

    .progress_msg(progress, "verbose", "    ", nrow(d), " records after quality filters")
    if (nrow(d) == 0) {
      return(tibble::tibble(species = clean_name, rows = 0L,
                            file_path = NA_character_, status = "empty_after_clean"))
    }

    ### H3 thinning: one random record per hex cell; tot_obs = distinct records in it.
    d_thin <- d |>
      dplyr::mutate(
        h3_index = suppressMessages(h3jsr::point_to_cell(
          data.frame(lng = decimalLongitude, lat = decimalLatitude), res = cfg$hex_resolution))) |>
      dplyr::mutate(tot_obs = dplyr::n_distinct(record_id), .by = h3_index) |>
      dplyr::slice_sample(prop = 1) |>
      dplyr::filter(!duplicated(h3_index)) |>
      dplyr::relocate(decimalLongitude, decimalLatitude, .after = dplyr::last_col())

    .progress_msg(progress, "verbose", "    ", nrow(d_thin), " records after H3 thinning")
    out_path <- file.path(ensure_dir(output_dir), paste0(clean_name, "_clean.parquet"))
    arrow::write_parquet(d_thin, out_path)
    tibble::tibble(species = clean_name, rows = nrow(d_thin),
                   file_path = as.character(fs::path_rel(out_path, start = here())),
                   status = "written")
  }, error = function(e) {
    tibble::tibble(species = clean_name, rows = NA_integer_,
                   file_path = NA_character_,
                   status = paste0("error: ", conditionMessage(e)))
  })
  res
}

### --- Cache -------------------------------------------------------------------
clean_manifest_path <- function(taxon, year, cfg = config) {
  file.path(taxon_output_dir(taxon, cfg), "logs",
            glue("clean_occurrence_manifest_{taxon}_{year}.json"))
}

### Hash of the cleaning code (deparsed, so comments/whitespace don't matter) and of
### every setting that affects the output.
clean_settings_hash <- function(cfg = config) {
  fns <- c(".normalise_gbif_cols", ".ref_buffers", "cc_near_ref", "clean_one_file")
  rlang::hash(list(
    code       = lapply(fns, function(f) deparse(get(f, envir = globalenv()))),
    cleaning   = cfg$cleaning,
    resolution = cfg$hex_resolution,
    crs        = cfg$crs$working
  ))
}

.read_clean_manifest <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

### Clean one taxon-year: every parsed occ file in that folder, skipping files whose
### input md5 and settings match the manifest and whose output still exists.
clean_taxon_year <- function(taxon, year, cfg = config, force = FALSE,
                             progress = c("quiet", "terse", "verbose")) {
  progress <- match.arg(progress)
  in_dir <- taxon_year_occ_dir(taxon, year, cfg)
  files <- list.files(in_dir, pattern = "_occs\\.parquet$", full.names = TRUE)
  if (length(files) == 0) return(tibble::tibble())

  settings <- clean_settings_hash(cfg)
  mpath <- clean_manifest_path(taxon, year, cfg)
  old <- if (force) NULL else .read_clean_manifest(mpath)
  old_files <- if (!is.null(old) && identical(old$settings, settings)) old$files else NULL

  md5 <- unname(tools::md5sum(files))
  .progress_msg(progress, "terse", glue("cleaning {taxon} for {year} ({length(files)} files)"))
  results <- purrr::pmap_dfr(list(files, md5, seq_along(files)), function(f, h, i) {
    nm <- sub("_occs\\.parquet$", "", basename(f))
    tag <- glue("  [{i}/{length(files)}] {nm}")
    prev <- if (!is.null(old_files)) old_files[old_files$input_md5 == h, ] else NULL
    if (!is.null(prev) && nrow(prev) > 0) {
      prev <- prev[1, ]
      ok <- is.na(prev$file_path) || file.exists(here::here(prev$file_path))
      if (prev$status %in% c("written", "empty_after_clean") && ok) {
        .progress_msg(progress, "verbose", tag, ": cached")
        return(tibble::tibble(species = prev$species, rows = as.integer(prev$rows),
                              file_path = prev$file_path, status = prev$status, cached = TRUE))
      }
    }
    .progress_msg(progress, "verbose", tag, ": cleaning")
    clean_one_file(f, output_dir = in_dir, cfg = cfg, progress = progress) |>
      dplyr::mutate(cached = FALSE)
  })

  ### Manifest last; only successful outcomes are cached (errors always retry).
  keep <- results$status %in% c("written", "empty_after_clean")
  ensure_parent_dir(mpath)
  jsonlite::write_json(
    list(taxon = taxon, year = year, settings = settings,
         files = dplyr::mutate(results[keep, c("species", "rows", "file_path", "status")],
                               input_md5 = md5[keep])),
    mpath, auto_unbox = TRUE, pretty = TRUE, na = "null"
  )

  results |> dplyr::mutate(taxon = taxon, year = year, .before = 1)
}

run_clean_occurrence <- function(taxa = config$taxa, years = config$years, cfg = config,
                                 force = FALSE,
                                 progress = c("terse", "verbose", "quiet")) {
  progress <- match.arg(progress)
  grid <- tidyr::expand_grid(taxon = taxa, year = years)
  results <- purrr::pmap_dfr(list(grid$taxon, grid$year),
                             function(taxon, year) clean_taxon_year(taxon, year, cfg, force, progress))
  log_path <- file.path(cfg$paths$output, "indicators", "species_populations",
                        "logs", "clean_occurrence.csv")
  ensure_parent_dir(log_path)
  readr::write_csv(results, log_path)
  results
}
