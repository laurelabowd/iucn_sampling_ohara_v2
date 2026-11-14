### =============================================================================
### R/03_process/build_hex_grid.R  —  region H3 grid + conservation-area cells
### =============================================================================
###
### Produces the two shared spatial layers the rest of the pipeline joins on:
###   hex_grid_path(config)        — hex_grid_res<N>.gpkg: every H3 cell covering the region
###   cons_areas_hex_path(config)  — cons_areas_hex_res<N>.gpkg: the subset of cells inside conservation areas,
###                                  with managing-agency / GAP attributes.
###
### Reworked from get_hex_grid.R: region and resolution come from config; the
### conservation-layer field names are configurable; the bird eBird-raster zonal
### block (a range source, not a grid step) is moved to range_from_raster.R; and
### the mangled double-assignment bug at get_hex_grid.R:144-147 is dropped.
###
### Requires config + R/00_utils.R.  Usage:  build_hex_grid()

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(h3jsr)
  library(rnaturalearth)
  library(here)
})

### Load the study-region boundary as an sf polygon in the working CRS.
get_region_boundary <- function(cfg = config) {
  ### options for boundary_source: 
  ### * "global" for entire globe
  ### * bbox in degrees (4-element vector, optional named xmin, xmax, ymin, ymax)
  ### * "rnaturalearth" with region name to pull polygon
  ### * a path to a spatial file readable by st_read()
  src <- cfg$region$boundary_source
  
  ### set up boundary_source = "global" option
  if (identical(src, "global")) {
    ### Return world bounds (entire globe)
    boundary <- sf::st_sf(
      name = "World",
      geometry = sf::st_sfc(
        sf::st_polygon(list(rbind(
          c(-180, -90), c(180, -90), c(180, 90), c(-180, 90), c(-180, -90)
        ))),
        crs = 4326
      )
    )
    return(boundary)
  }

  ### set up for boundary_source = "bbox" in deg xmin/xmax/ymin/ymax
  if (identical(src, "bbox")) {
    ### Create polygon from bounding box
    bbox <- cfg$region$bbox

    ### if bbox vec is named, use names in whatever order given; 
    ### otherwise set names in assumed order
    bbox_names <- c('xmin', 'xmax', 'ymin', 'ymax')
    if(length(bbox) != 4) stop('bbox must include 4 values: xmin, xmax, ymin, ymax in degrees')
    if(any(!names(bbox) %in% bbox_names)) stop('bbox must include 4 values: xmin, xmax, ymin, ymax in degrees')
    if(is.null(names(bbox))) names(bbox) <- bbox_names
    polygon <- sf::st_polygon(list(rbind(
      c(bbox['xmin'], bbox['ymin']),
      c(bbox['xmin'], bbox['ymax']),
      c(bbox['xmax'], bbox['ymax']),
      c(bbox['xmax'], bbox['ymin']),
      c(bbox['xmin'], bbox['ymin'])  ### close ring
    )))
    boundary <- sf::st_sf(name = cfg$region$admin_name, geometry = sf::st_sfc(polygon, crs = 4326))
    return(boundary)
  }

  ### set up boundary_source = "rnaturalearth" and filter by state
  if (identical(src, "rnaturalearth")) {
    b <- rnaturalearth::ne_states(country = "United States of America") |>
      ### update: name %in% cfg$region$admin_name so multiple regions
      ### can be added?
      dplyr::filter(name %in% cfg$region$admin_name) |>
      dplyr::select(name)
  } else {
    b <- sf::st_read(src, quiet = TRUE)
  }
  b_out <- sf::st_transform(b, cfg$crs$working)
  return(b_out)
}

### Load and filter the conservation-areas layer to strictly-protected terrestrial
### units, dissolving to one geometry per management unit.
get_conservation_areas <- function(cfg = config) {
  ca <- cfg$conservation_areas
  sfp <- (if (is.null(ca$layer)) sf::st_read(ca$path, quiet = TRUE)
          else sf::st_read(ca$path, layer = ca$layer, quiet = TRUE)) |>
    sf::st_zm(drop = TRUE)

  needed <- c(ca$gap_field, ca$realm_field, ca$unit_field, ca$agency_field, ca$agency_level_field)
  missing_cols <- setdiff(needed, names(sfp))
  if (length(missing_cols)) {
    stop("Conservation layer is missing field(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  keep <- sfp[[ca$gap_field]] %in% ca$gap_values
  if (!is.null(ca$realm_field)) keep <- keep & sfp[[ca$realm_field]] %in% ca$realm_value
  sfp <- sfp[keep, ]
  if (nrow(sfp) == 0) stop("No conservation areas left after gap/realm filtering.", call. = FALSE)

  ### No dissolve: a union needs valid geometry, and per-cell attributes are
  ### de-duplicated downstream anyway. s2 requires lon/lat.
  sfp |>
    dplyr::select(dplyr::all_of(
      c(ca$unit_field, ca$gap_field, ca$agency_field, ca$agency_level_field))) |>
    sf::st_transform(4326)
}

### --- Build manifest ----------------------------------------------------------
### Records what the grid layers were built from so an unchanged request is not
### rebuilt (~5 min at level 4). Written next to the layers, per H3 level.
hex_manifest_path <- function(cfg = config) {
  hex_layer_path(file.path(dirname(cfg$paths$hex_grid), "hex_grid_manifest.json"), cfg)
}

### Hash of everything that determines the layers: region, resolution, the
### conservation-area settings, and the source file's names, sizes and mtimes
### (a .gdb is a directory, so every file inside is included).
hex_grid_request_hash <- function(cfg = config) {
  src <- cfg$conservation_areas$path
  files <- if (dir.exists(src)) list.files(src, recursive = TRUE, full.names = TRUE) else src
  fi <- file.info(files)
  rlang::hash(list(
    region     = cfg$region,
    resolution = cfg$hex_resolution,
    cons       = cfg$conservation_areas,
    source     = data.frame(file = basename(files), size = fi$size,
                            mtime = as.numeric(fi$mtime))
  ))
}

### Returns NULL if existing layers are reusable, else a short reason string.
hex_grid_stale_reason <- function(cfg = config) {
  if (!file.exists(hex_grid_path(cfg)) || !file.exists(cons_areas_hex_path(cfg))) {
    return("layer file missing")
  }
  p <- hex_manifest_path(cfg)
  if (!file.exists(p)) return("no manifest")
  m <- tryCatch(jsonlite::read_json(p, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(m)) return("unreadable manifest")
  if (!identical(m$request_hash, hex_grid_request_hash(cfg))) return("inputs changed")
  NULL
}

### Build (or reuse) the region hex grid and conservation-area cells.
### force = TRUE always rebuilds.
build_hex_grid <- function(cfg = config, force = FALSE) {
  if (!isTRUE(force)) {
    reason <- hex_grid_stale_reason(cfg)
    if (is.null(reason)) {
      message(glue::glue("[hex_grid] level {cfg$hex_resolution}: using CACHED layers ",
                         "(inputs unchanged; force = TRUE to rebuild)."))
      return(invisible(list(
        hex_grid   = sf::st_read(hex_grid_path(cfg), quiet = TRUE),
        cons_cells = sf::st_read(cons_areas_hex_path(cfg), quiet = TRUE)
      )))
    }
    message(glue::glue("[hex_grid] level {cfg$hex_resolution}: building ({reason})."))
  }

  region <- get_region_boundary(cfg)

  ### Region H3 grid.
  ### polygon_to_cells() cannot fill a globe-spanning polygon (returns NA), so
  ### global regions are built from res-0 cells and their children instead.
  bb <- sf::st_bbox(sf::st_transform(region, 4326))
  is_global <- identical(cfg$region$boundary_source, "global") ||
    (bb[["xmin"]] <= -179.9 && bb[["xmax"]] >= 179.9 &&
       bb[["ymin"]] <= -89.9 && bb[["ymax"]] >= 89.9)

  ids <- if (is_global) {
    h3jsr::get_children(h3jsr::get_res0(), res = cfg$hex_resolution, simple = TRUE)
  } else {
    h3jsr::polygon_to_cells(region, res = cfg$hex_resolution, simple = TRUE)
  }
  ids <- unique(unlist(ids, use.names = FALSE))
  if (length(ids) == 0 || anyNA(ids)) {
    stop("polygon_to_cells() returned no valid H3 cells; the region may be too ",
         "large or not in lon/lat (EPSG:4326).", call. = FALSE)
  }
  hex_grid <- h3jsr::cell_to_polygon(ids, simple = FALSE) |>
    dplyr::rename(h3_index = h3_address)
  ensure_parent_dir(hex_grid_path(cfg))
  sf::st_write(hex_grid, hex_grid_path(cfg), delete_dsn = TRUE, quiet = TRUE)

  ### Conservation-area cells (attributes summarised per cell).
  ca_cfg <- cfg$conservation_areas
  cons <- get_conservation_areas(cfg)

  ### Spherical intersects via s2. WDPA polygons can be s2-invalid even after
  ### st_make_valid() (GEOS validity is planar), so rebuild them with split
  ### crossing edges. Planar joins give false hits at the antimeridian.
  cons_s2 <- s2::as_s2_geography(sf::st_geometry(cons), oriented = FALSE, check = FALSE) |>
    s2::s2_rebuild(options = s2::s2_options(
      edge_type = "undirected", validate = FALSE, split_crossing_edges = TRUE))
  hex_s2 <- s2::as_s2_geography(sf::st_geometry(hex_grid))
  hits <- s2::s2_intersects_matrix(hex_s2, cons_s2)

  pairs <- data.frame(
    h3_index = rep(hex_grid$h3_index, lengths(hits)),
    cons_row = unlist(hits, use.names = FALSE)
  )
  cons_attr <- sf::st_drop_geometry(cons)[pairs$cons_row, , drop = FALSE]

  cell_attrs <- dplyr::bind_cols(pairs["h3_index"], cons_attr) |>
    dplyr::group_by(h3_index) |>
    dplyr::summarize(
      unit_name  = paste(unique(na.omit(.data[[ca_cfg$unit_field]])), collapse = "; "),
      mng_agency = paste(unique(na.omit(.data[[ca_cfg$agency_field]])), collapse = "; "),
      gap_status = paste(sort(unique(na.omit(.data[[ca_cfg$gap_field]]))), collapse = ", "),
      .groups = "drop"
    )
  cons_cells <- dplyr::inner_join(hex_grid, cell_attrs, by = "h3_index")
  sf::st_write(cons_cells, cons_areas_hex_path(cfg), delete_dsn = TRUE, quiet = TRUE)

  ### Manifest is written last, so an interrupted build is never treated as cached.
  jsonlite::write_json(
    list(request_hash = hex_grid_request_hash(cfg), hex_resolution = cfg$hex_resolution,
         n_hex_cells = nrow(hex_grid), n_cons_cells = nrow(cons_cells),
         built_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S")),
    hex_manifest_path(cfg), auto_unbox = TRUE, pretty = TRUE
  )

  message(glue::glue("[hex_grid] {nrow(hex_grid)} region cells, ",
                     "{nrow(cons_cells)} conservation cells written."))
  invisible(list(hex_grid = hex_grid, cons_cells = cons_cells))
}
