### =============================================================================
### R/03_process/04c_range_from_polygons_h3o_fn.R  —  range step using h3o (experimental)
### =============================================================================
###
### Parallel implementation of range_from_polygons using h3o instead of h3jsr.
### h3o is an R wrapper around the Rust h3 library and may be faster for large
### polygon-to-cells conversions.
###
### Function signature mirrors 04a_range_from_polygons_fn.R; designed for
### drop-in replacement and benchmarking against h3jsr.
###
### Requires config + taxa_config + R/00_utils.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(glue)
  library(sf)
  library(h3o)
  library(arrow)
  library(here)
})

### Map each range polygon -> its canonical species name via the crosswalk.
### (identical to h3jsr version)
.polygon_species_map <- function(taxon, tcfg, cfg = config) {
  cw <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE)
  if (!"range_polygon_layer" %in% names(cw)) cw$range_polygon_layer <- NA_character_
  cw <- cw |>
    dplyr::filter(taxon == !!taxon, !is.na(range_polygon_name), !is.na(name_gbif)) |>
    dplyr::transmute(range_polygon_name, species = name_gbif,
                     polygon_path = range_polygon_path,
                     polygon_layer = as.character(range_polygon_layer)) |>
    tidyr::separate_longer_delim(c(range_polygon_name, polygon_path, polygon_layer),
                                 delim = "|")
  if (nrow(cw) == 0) {
    warning(glue("[{taxon}] no polygon-linked species in crosswalk."))
  }
  cw
}

### Polygon area in km2 (NA if the geometry can't be measured).
### (identical to h3jsr version)
.area_km2 <- function(sf_obj) {
  tryCatch(as.numeric(sum(sf::st_area(sf::st_make_valid(sf_obj)))) / 1e6,
           error = function(e) NA_real_)
}

### H3 cells for a polygon set using h3o.
### h3o::polygon_to_cells is simpler and doesn't fail on antimeridian-crossing
### polygons the way h3jsr::polygon_to_cells does. If it still fails, fall back
### to strip-cutting logic.
.polygon_cells_h3o <- function(sp_sf, res) {

  ### attempt to convert polygon to cells using h3o
  ### (h3o has no polygon_to_cells(); sfc_to_cells + centroid matches h3jsr)
  to_cells <- function(x) {
    h3o::sfc_to_cells(sf::st_geometry(x), resolution = res, containment = "centroid") |>
      h3o::flatten_h3() |>
      as.character()
  }
  direct <- tryCatch(to_cells(sp_sf), error = function(e) NULL)
  ### if !is.null, conversion was successful - return the result!
  if (!is.null(direct)) return(direct)

  ### if not successful, fall back to strip-cutting (same as h3jsr version)
  ### This handles edge cases like invalid geometries.
  old_s2 <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  geom <- sf::st_make_valid(sf::st_geometry(sp_sf))
  sf::st_crs(geom) <- 4326
  strips <- lapply(seq(-180, 120, by = 60), function(x0) {
    sf::st_sfc(sf::st_polygon(list(rbind(c(x0, -90), c(x0 + 60, -90), c(x0 + 60, 90),
                                         c(x0, 90), c(x0, -90)))), crs = 4326)
  })
  unique(unlist(lapply(strips, function(s) {
    piece <- suppressMessages(suppressWarnings(sf::st_intersection(geom, s)))
    piece <- piece[!sf::st_is_empty(piece)]
    if (length(piece) == 0) return(NULL)
    suppressMessages(to_cells(sf::st_sf(geometry = piece)))
  })))
}

### Read one polygon layer once and return H3 cells for each requested species.
### (Identical to h3jsr version except calls .polygon_cells_h3o instead of .polygon_cells)
.cells_for_layer_h3o <- function(path, layer, species_names, name_field, cfg = config,
                                 progress = "terse", taxon = "") {

  lyr_lab <- if (is.na(layer)) "default layer" else glue("layer '{layer}'")
  .progress_msg(progress, "verbose",
                glue("  [{taxon}] reading {basename(path)} ({lyr_lab}, ",
                     "{nrow(species_names)} species to extract)"))
  range_sf <- if (is.na(layer)) sf::st_read(path, quiet = TRUE)
              else sf::st_read(path, layer = layer, quiet = TRUE)
  if (!name_field %in% names(range_sf)) {
    stop(glue("no field '{name_field}' in {basename(path)}"))
  }
  .progress_msg(progress, "verbose",
                glue("  [{taxon}] {nrow(range_sf)} features read; transforming to working CRS"))
  range_sf <- sf::st_transform(range_sf, cfg$crs$working)

  purrr::pmap_dfr(species_names, function(range_polygon_name, species) {
    sp_sf <- range_sf[as.character(range_sf[[name_field]]) %in% range_polygon_name, ]
    
    if (nrow(sp_sf) == 0) return(NULL)

    n_feat <- nrow(sp_sf)
    detail <- if (progress == "verbose") {
      glue(" ({n_feat} feature{if (n_feat > 1) 's' else ''}, {format(round(.area_km2(sp_sf)), big.mark = ',')} km2)")
    } else if (n_feat > 1) {
      glue(" ({n_feat} features combined)")
    } else ""
    .progress_msg(progress, "verbose", glue("    {species}{detail}"))

    t0 <- Sys.time()
    cells <- .polygon_cells_h3o(sp_sf, cfg$hex_resolution) |> unlist()
    cells <- as.character(unique(cells))
    .progress_msg(progress, "verbose",
                  glue("      {length(cells)} cells in ",
                       "{round(as.numeric(difftime(Sys.time(), t0, units = 'secs')), 1)}s"))
    tibble::tibble(species = species, h3_index = cells,
                   n_features = n_feat, source = paste(path, layer))
  })
}

### Build the wide range-cell table for one taxon using h3o (drop-in replacement for
### build_range_from_polygons with h3jsr).
build_range_from_polygons_h3o <- function(taxon, cfg = config, taxa_cfg = taxa_config,
                                           samp_cfg = NULL,
                                           progress = c("terse", "verbose", "quiet")) {
  progress <- match.arg(progress)
  tcfg <- taxa_cfg[[taxon]]
  if (!is_polygon_source(tcfg$range_source)) {
    .progress_msg(progress, "terse",
                  glue("[{taxon}] range_source is '{tcfg$range_source}', not a polygon source; skipped."))
    return(invisible(NULL))
  }
  if (is.null(samp_cfg)) samp_cfg <- make_sampling_config(cfg$focal_year, cfg = cfg)[[taxon]]

  hex_grid <- sf::st_read(hex_grid_path(cfg), quiet = TRUE) |> sf::st_drop_geometry()
  cons_ids <- sf::st_read(cons_areas_hex_path(cfg), quiet = TRUE)$h3_index

  spp_map <- .polygon_species_map(taxon, tcfg, cfg)

  n_sources <- dplyr::distinct(spp_map, polygon_path, polygon_layer) |> nrow()
  n_files   <- dplyr::n_distinct(spp_map$polygon_path)
  .progress_msg(progress, "terse",
                glue("[{taxon}] building ranges for {dplyr::n_distinct(spp_map$species)} species ",
                     "from {n_files} file{if (n_files != 1) 's' else ''}",
                     "{if (n_sources > n_files) glue(' / {n_sources} layers') else ''}"))

  multi_src <- spp_map |>
    dplyr::summarise(n_src = dplyr::n(), .by = species) |>
    dplyr::filter(n_src > 1)
  if (nrow(multi_src) > 0) {
    .progress_msg(progress, "terse",
                  glue("[{taxon}] {nrow(multi_src)} species have polygons in >1 layer/file ",
                       "(combined): {paste(utils::head(multi_src$species, 5), collapse = ', ')}",
                       "{if (nrow(multi_src) > 5) ', ...' else ''}"))
  }

  long <- spp_map |>
    dplyr::group_split(polygon_path, polygon_layer, .keep = TRUE) |>
    purrr::map_dfr(function(g) {
      path <- g$polygon_path[1]; layer <- g$polygon_layer[1]
      tryCatch(
        .cells_for_layer_h3o(path, layer, g[c("range_polygon_name", "species")],
                             tcfg$polygon_name_field, cfg,
                             progress = progress, taxon = taxon),
        error = function(e) {
          warning(glue("failed {basename(path)} [{layer}]: {conditionMessage(e)}"))
          NULL
        }
      )
    })

  if (nrow(long) > 0) {
    multi <- long |>
      dplyr::distinct(species, source, n_features) |>
      dplyr::summarise(n_features = sum(n_features), .by = species) |>
      dplyr::filter(n_features > 1)
    if (nrow(multi) > 0) {
      .progress_msg(progress, "terse",
                    glue("[{taxon}] {nrow(multi)} species combine >1 polygon feature: ",
                         "{paste(utils::head(multi$species, 5), collapse = ', ')}",
                         "{if (nrow(multi) > 5) ', ...' else ''}",
                         " (see spp_multi_feature_ranges.csv)"))
      ensure_dir(taxon_output_dir(taxon, cfg))
      readr::write_csv(multi, file.path(taxon_output_dir(taxon, cfg),
                                        "spp_multi_feature_ranges_h3o.csv"))
    }
    long <- dplyr::select(long, species, h3_index)
  }

  if (nrow(long) == 0) {
    stop(glue("[{taxon}] no range cells were built for any species; see warnings above."),
         call. = FALSE)
  }
  wide <- long |>
    dplyr::mutate(presence = 1L) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(names_from = species, values_from = presence, values_fill = 0L)

  range_cells <- hex_grid |>
    dplyr::left_join(wide, by = "h3_index") |>
    dplyr::mutate(dplyr::across(-h3_index, ~tidyr::replace_na(., 0L)))

  invisible(range_cells)
}
