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
  library(h3o)
  library(arrow)
  library(here)
})

### Settings
### Northern latitude cap (degrees) applied to range polygons before H3 conversion.
h3_max_lat <- 88

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

### Polygon area in km2 (NA if the geometry can't be measured). Display only (verbose
### progress). Planar area in the equal-area CRS: ~1 s for 3M vertices and within 0.02% of
### geodesic area on Delphinus/Phocoena, where s2 st_area errored on make_valid output
### ("Loop ... crosses loop") and geodesic area with s2 off took ~150 s.
.area_km2 <- function(sf_obj, crs = config$crs$equal_area) {
  tryCatch(as.numeric(sum(sf::st_area(sf::st_transform(sf::st_geometry(sf_obj), crs)))) / 1e6,
           error = function(e) NA_real_)
}

### H3 cells for a polygon set, using h3o (Rust; ~15x faster than h3jsr on abalones,
### identical cells - see logs/DECISIONS.md 2026-11-18). h3o has no polygon_to_cells();
### sfc_to_cells() is its equivalent. containment = "centroid" reproduces h3jsr, which
### keeps cells whose centre falls inside the polygon. All features passed in are
### flattened into one cell vector, so multi-feature species become one range.
### Circumpolar polygons (e.g. Ursus maritimus) hang h3o instead of erroring: H3 treats
### edges as great circles, so a flat edge along the top latitude spanning 360 degrees
### of longitude collapses to a sliver around the pole. A hang never reaches a
### tryCatch, so the polygons are first cropped to latitudes <= `max_lat` (the polar
### cap is not needed; see logs/DECISIONS.md), and any feature wider than 180 degrees
### of longitude goes straight to 60-degree longitude strips with the cells united.
### Narrower features try the direct conversion, with the strips as the error fallback.
### Planar clipping, hence s2 is switched off.
.polygon_cells <- function(sp_sf, res, max_lat = h3_max_lat) {

  to_cells <- function(x) {
    h3o::sfc_to_cells(sf::st_geometry(x), resolution = res, containment = "centroid") |>
      h3o::flatten_h3() |>
      as.character()
  }

  ### st_intersection can return stray lines/points (where a polygon only touches the
  ### clip edge) or mixed geometry types, which h3o rejects: keep polygon parts only.
  polys_only <- function(x) {
    x <- suppressWarnings(sf::st_collection_extract(x, "POLYGON", warn = FALSE))
    x <- x[!sf::st_is_empty(x)]
    x <- sf::st_cast(x, "MULTIPOLYGON")
    ### st_collection_extract leaves nested ALTREP wrap_list objects that h3o's Rust code
    ### cannot take a writable pointer to ("cannot take a writable DATAPTR"); a WKB round
    ### trip rebuilds plain lists.
    sf::st_as_sfc(sf::st_as_binary(x), crs = 4326)
  }

  old_s2 <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  geom <- sf::st_geometry(sp_sf)
  sf::st_crs(geom) <- 4326

  ### crop the north polar cap (northern edge only; the south is left as is).
  ### st_make_valid is slow on big multipolygons (~2.5 min for 3M vertices) and often
  ### unnecessary, so it is only used if the planar crop errors on invalid geometry.
  crop_box <- sf::st_sfc(sf::st_polygon(list(rbind(c(-180, -90), c(180, -90), c(180, max_lat),
                                                   c(-180, max_lat), c(-180, -90)))), crs = 4326)
  crop <- function(g) suppressMessages(suppressWarnings(sf::st_intersection(g, crop_box)))
  geom <- tryCatch(crop(geom), error = function(e) crop(sf::st_make_valid(geom))) |>
    polys_only()
  if (length(geom) == 0) return(character(0))

  ### Explode to single polygons: a multipolygon's bbox can span 360 degrees only because
  ### of scattered islands, and each polygon is what h3o has to fill. Only polygons wider
  ### than 180 degrees of longitude (circumpolar) need strips.
  geom <- sf::st_cast(geom, "POLYGON")
  lon_span <- vapply(geom, function(p) diff(range(p[[1]][, 1])), numeric(1))
  narrow <- lon_span <= 180

  cells <- character(0)
  if (any(narrow)) {
    ### narrow polygons: attempt direct conversion; on failure send everything to strips
    direct <- tryCatch(to_cells(polys_only(geom[narrow])), error = function(e) NULL)
    if (is.null(direct)) narrow[] <- FALSE else cells <- direct
  }
  if (all(narrow)) return(unique(cells))

  ### wide polygons (or direct conversion failed): cut to strips and convert each
  wide <- geom[!narrow]
  strips <- lapply(seq(-180, 120, by = 60), function(x0) {
    sf::st_sfc(sf::st_polygon(list(rbind(c(x0, -90), c(x0 + 60, -90), c(x0 + 60, 90),
                                         c(x0, 90), c(x0, -90)))), crs = 4326)
  })
  unique(c(cells, unlist(lapply(strips, function(s) {
    piece <- suppressMessages(suppressWarnings(sf::st_intersection(wide, s))) |>
      polys_only()
    if (length(piece) == 0) return(NULL)
    to_cells(sf::st_sf(geometry = piece))
  }))))
}

### Read one polygon layer once and return H3 cells for each requested species.
### A layer may hold many species (rows matched on `name_field`) or just one.
### `species_names` maps polygon names to canonical species: a data frame with
### range_polygon_name and species. Returns species, h3_index, n_features, where
### n_features is the number of polygon features combined for that species here.
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

    ### Convert the polygons to H3 cells, timing it
    t0 <- Sys.time() ### (only used if "verbose")
    cells <- .polygon_cells(sp_sf, cfg$hex_resolution)
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
