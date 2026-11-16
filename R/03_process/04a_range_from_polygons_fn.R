### =============================================================================
### R/03_process/04a_range_from_polygons_fn.R  —  CANONICAL range step (functions)
### =============================================================================
###
### Turn range polygons (multi-species layers or one file per species) into a wide table of H3 range cells:
### one row per hex cell, one 0/1 column per species (keyed by the GBIF canonical
### name), plus per-species totals overall and within conservation areas.
###
### Reworked from get_hex_grid_from_shape.R, which overwrote its polygon list
### three times (so only mammals ever ran) and hardcoded the mammal output path.
### Here it is a function of `taxon`, looped by the _run driver.
###
### The species identity comes from the crosswalk, not the raw polygon: each
### polygon's name field (e.g. sci_name) was matched to name_gbif in Stage 0, so range
### cells and occurrence/richness data share ONE canonical key and join cleanly.
###
### Requires config + taxa_config + R/00_utils.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(glue)
  library(sf)
  library(h3jsr)
  library(arrow)
  library(here)
})

### Map each range polygon -> its canonical species name via the crosswalk.
### One row per (species, source). Returns a tibble: polygon_path, polygon_layer, range_polygon_name,
### species (name_gbif). polygon_layer is NA for crosswalks built before layers
### were recorded (read as the file's default layer).
.polygon_species_map <- function(taxon, tcfg, cfg = config) {
  cw <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE)
  if (!"range_polygon_layer" %in% names(cw)) cw$range_polygon_layer <- NA_character_
  cw <- cw |>
    dplyr::filter(taxon == !!taxon, !is.na(range_polygon_name), !is.na(name_gbif)) |>
    dplyr::transmute(range_polygon_name, species = name_gbif,
                     polygon_path = range_polygon_path,
                     polygon_layer = as.character(range_polygon_layer)) |>
    ### Species with several sources hold "|"-delimited parallel values.
    tidyr::separate_longer_delim(c(range_polygon_name, polygon_path, polygon_layer),
                                 delim = "|")
  if (nrow(cw) == 0) {
    warning(glue("[{taxon}] no polygon-linked species in crosswalk."))
  }
  cw
}

### Polygon area in km2 (NA if the geometry can't be measured).
.area_km2 <- function(sf_obj) {
  tryCatch(as.numeric(sum(sf::st_area(sf::st_make_valid(sf_obj)))) / 1e6,
           error = function(e) NA_real_)
}

### H3 cells for a polygon set. polygon_to_cells() fails (code 1) on polygons that
### span most of the globe in longitude (e.g. circumpolar Southern Ocean ranges with
### edges on the antimeridian), so on failure the polygons are cut into 60-degree
### longitude strips and the cells united. Planar clipping, hence s2 is switched off.
.polygon_cells <- function(sp_sf, res) {

  ### attempt to convert polygon to cells
  direct <- tryCatch(h3jsr::polygon_to_cells(sp_sf, res = res), error = function(e) NULL)
  ### if !is.null, conversion was successful - return the result!
  if (!is.null(direct)) return(direct)

  ### if not successful, could be invalid geom, or could be too
  ### wide of polygon: convert to planar geom, make_valid, cut to strips,
  ### then convert each strip to polygons and rejoin
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
    suppressMessages(h3jsr::polygon_to_cells(sf::st_sf(geometry = piece), res = res))
  })))
}

### Read one polygon layer once and return H3 cells for each requested species.
### A layer may hold many species (rows matched on `name_field`) or just one.
### `species_names` maps polygon names to canonical species: a data frame with
### range_polygon_name and species. Returns species, h3_index, n_features, where
### n_features is the number of polygon features combined for that species here.
### h3jsr's own "this might take a while" message is replaced by a per-species
### progress line (terse) so it is attributed to the species being processed.
.cells_for_layer <- function(path, layer, species_names, name_field, cfg = config,
                             progress = "terse", taxon = "") {

  ### Announce which file/layer is being read and how many species we want from it
  lyr_lab <- if (is.na(layer)) "default layer" else glue("layer '{layer}'")
  .progress_msg(progress, "verbose",
                glue("  [{taxon}] reading {basename(path)} ({lyr_lab}, ",
                     "{nrow(species_names)} species to extract)"))
  ### Read the layer once (default layer if none was recorded in the crosswalk)
  range_sf <- if (is.na(layer)) sf::st_read(path, quiet = TRUE)
              else sf::st_read(path, layer = layer, quiet = TRUE)
  ### Fail early (caught by the caller as a warning) if the species-name column is absent
  if (!name_field %in% names(range_sf)) {
    stop(glue("no field '{name_field}' in {basename(path)}"))
  }
  ### Reproject once so every species shares the working CRS
  .progress_msg(progress, "verbose",
                glue("  [{taxon}] {nrow(range_sf)} features read; transforming to working CRS"))
  range_sf <- sf::st_transform(range_sf, cfg$crs$working)

  ### Process each requested species in turn and stack the results into one table
  purrr::pmap_dfr(species_names, function(range_polygon_name, species) {
    ### Pull out all polygon features whose name matches this species
    sp_sf <- range_sf[as.character(range_sf[[name_field]]) %in% range_polygon_name, ]
    
    ### if no rows in this layer, spp not present in this layer; return empty value
    if (nrow(sp_sf) == 0) return(NULL)

    ### Build a short description for progress messages: feature count (+ area if verbose)
    n_feat <- nrow(sp_sf)
    detail <- if (progress == "verbose") {
      glue(" ({n_feat} feature{if (n_feat > 1) 's' else ''}, {format(round(.area_km2(sp_sf)), big.mark = ',')} km2)")
    } else if (n_feat > 1) {
      glue(" ({n_feat} features combined)")
    } else ""
    .progress_msg(progress, "verbose", glue("    {species}{detail}"))

    ### Convert the polygons to H3 cells, timing it. Large ranges are slow, so the
    ### h3jsr "might take a while" message is replaced with a species-specific one.
    t0 <- Sys.time() ### (only used if "verbose")
    cells <- withCallingHandlers(
      .polygon_cells(sp_sf, cfg$hex_resolution) |> unlist(),
      message = function(m) {
        if (grepl("might take a while", conditionMessage(m), fixed = TRUE)) {
          .progress_msg(progress, "terse",
                        glue("  [{taxon}] {species}{detail}: large range, this might take a while..."))
          invokeRestart("muffleMessage")
        }
      }
    )
    ### Drop duplicate cells (overlapping features) and make sure they are plain strings
    cells <- as.character(unique(cells))
    .progress_msg(progress, "verbose",
                  glue("      {length(cells)} cells in ",
                       "{round(as.numeric(difftime(Sys.time(), t0, units = 'secs')), 1)}s"))
    # One row per cell for this species; n_features and source support later diagnostics
    tibble::tibble(species = species, h3_index = cells,
                   n_features = n_feat, source = paste(path, layer))
  })
}

### Build the wide range-cell table for one taxon and write it + summaries.
build_range_from_polygons <- function(taxon, cfg = config, taxa_cfg = taxa_config,
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

  ### Full CA hex grid and conservation-area cells (from 03_build_hex_grid.R).
  hex_grid <- sf::st_read(hex_grid_path(cfg), quiet = TRUE) |> sf::st_drop_geometry()
  cons_ids <- sf::st_read(cons_areas_hex_path(cfg), quiet = TRUE)$h3_index

  spp_map <- .polygon_species_map(taxon, tcfg, cfg)

  n_sources <- dplyr::distinct(spp_map, polygon_path, polygon_layer) |> nrow()
  n_files   <- dplyr::n_distinct(spp_map$polygon_path)
  .progress_msg(progress, "terse",
                glue("[{taxon}] building ranges for {dplyr::n_distinct(spp_map$species)} species ",
                     "from {n_files} file{if (n_files != 1) 's' else ''}",
                     "{if (n_sources > n_files) glue(' / {n_sources} layers') else ''}"))

  ### Species whose polygons are spread over several layers/files are combined.
  multi_src <- spp_map |>
    dplyr::summarise(n_src = dplyr::n(), .by = species) |>
    dplyr::filter(n_src > 1)
  if (nrow(multi_src) > 0) {
    .progress_msg(progress, "terse",
                  glue("[{taxon}] {nrow(multi_src)} species have polygons in >1 layer/file ",
                       "(combined): {paste(utils::head(multi_src$species, 5), collapse = ', ')}",
                       "{if (nrow(multi_src) > 5) ', ...' else ''}"))
  }

  ### One read per (file, layer); species within it are split by name field.
  long <- spp_map |>
    dplyr::group_split(polygon_path, polygon_layer, .keep = TRUE) |>
    purrr::map_dfr(function(g) {
      path <- g$polygon_path[1]; layer <- g$polygon_layer[1]
      tryCatch(
        .cells_for_layer(path, layer, g[c("range_polygon_name", "species")],
                         tcfg$polygon_name_field, cfg,
                         progress = progress, taxon = taxon),
        error = function(e) {
          warning(glue("failed {basename(path)} [{layer}]: {conditionMessage(e)}"))
          NULL
        }
      )
    })

  ### Flag species built from more than one polygon feature (e.g. separate
  ### presence/seasonal polygons); all are combined into one range.
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
                                        "spp_multi_feature_ranges.csv"))
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

  ### Join onto the full grid so every cell is represented.
  range_cells <- hex_grid |>
    dplyr::left_join(wide, by = "h3_index") |>
    dplyr::mutate(dplyr::across(-h3_index, ~tidyr::replace_na(., 0L)))

  ensure_parent_dir(samp_cfg$range_cells_path)
  arrow::write_parquet(range_cells, samp_cfg$range_cells_path)

  ### Per-species totals: whole region and within conservation areas.
  tot <- range_cells |>
    dplyr::summarize(dplyr::across(-h3_index, ~sum(.x, na.rm = TRUE))) |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "species", values_to = "tot_range_cells")
  in_cons <- range_cells |>
    dplyr::filter(h3_index %in% cons_ids) |>
    dplyr::summarize(dplyr::across(-h3_index, ~sum(.x, na.rm = TRUE))) |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "species", values_to = "range_cells_in_cons_area")

  out_dir <- taxon_output_dir(taxon, cfg)
  ensure_dir(out_dir)
  readr::write_csv(tot, file.path(out_dir, "spp_cells_in_region.csv"))
  readr::write_csv(in_cons, file.path(out_dir, "spp_cells_in_cons_area.csv"))

  invisible(range_cells)
}
