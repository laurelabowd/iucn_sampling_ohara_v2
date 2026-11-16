### =============================================================================
### R/03_process/04_range_from_spatial_run.R  —  range-building driver
### =============================================================================
### Builds per-species hex range cells for every taxon, dispatching by
### range_source: "polygon" -> 04a_range_from_polygons_fn.R (aliases: "shapefile", "vector")
###               "raster"  -> 04b_range_from_raster_fn.R (eBird birds).
### `progress`: "terse" (default), "verbose" (per layer/species detail incl. area),
### or "quiet". Both write the same range_grid_cells.parquet shape, so downstream is uniform.
###
### Caching: per taxon, a manifest (range_manifest_<taxon>.json in the taxon's logs
### folder) records a hash of everything that determines the range outputs (crosswalk
### rows, range source files, hex grid layers, resolution, CRS, range code) plus the
### md5 of each output file (range_cells parquet and the two spp_cells_*.csv). A taxon
### is skipped when the hash matches and all outputs still exist with the recorded md5.
### The manifest is written last, so an interrupted build is never treated as cached.
### `force = TRUE` rebuilds everything.
###
### Requires config + taxa_config + sampling_config + R/00_utils.R.

source(here::here("R", "03_process", "04a_range_from_polygons_fn.R"))
source(here::here("R", "03_process", "04b_range_from_raster_fn.R"))

### --- Cache -------------------------------------------------------------------
range_manifest_path <- function(taxon, cfg = config) {
  file.path(taxon_output_dir(taxon, cfg), "logs", glue::glue("range_manifest_{taxon}.json"))
}

### Name, size and mtime of every file under `paths` (a .gdb is a directory).
.file_fingerprint <- function(paths) {
  paths <- unique(paths[!is.na(paths)])
  files <- unlist(lapply(paths, function(p) {
    if (dir.exists(p)) list.files(p, recursive = TRUE, full.names = TRUE) else p
  }))
  fi <- file.info(files)
  data.frame(file = files, size = fi$size, mtime = as.numeric(fi$mtime))[order(files), ]
}

.range_output_paths <- function(taxon, cfg = config) {
  out_dir <- taxon_output_dir(taxon, cfg)
  c(range_cells = make_sampling_config(cfg$focal_year, cfg = cfg)[[taxon]]$range_cells_path,
    cells_in_region = file.path(out_dir, "spp_cells_in_region.csv"),
    cells_in_cons_area = file.path(out_dir, "spp_cells_in_cons_area.csv"))
}

### Hash of everything that determines a taxon's range outputs.
range_settings_hash <- function(taxon, cfg = config, taxa_cfg = taxa_config) {
  tcfg <- taxa_cfg[[taxon]]
  src <- tcfg$range_source
  source_inputs <- if (is_polygon_source(src)) {
    spp_map <- suppressWarnings(.polygon_species_map(taxon, tcfg, cfg))
    list(map = spp_map, name_field = tcfg$polygon_name_field,
         files = .file_fingerprint(spp_map$polygon_path))
  } else if (identical(src, "raster")) {
    list(lut = .ebird_to_gbif(cfg),
         files = .file_fingerprint(tcfg$ebird_richness_raster))
  } else list(unknown = src)
  fns <- c("build_range_from_polygons", ".cells_for_layer", ".polygon_species_map",
           ".area_km2", ".polygon_cells", "build_range_from_raster", ".ebird_to_gbif")
  rlang::hash(list(
    taxon      = taxon,
    source     = src,
    inputs     = source_inputs,
    hex_grid   = unname(tools::md5sum(c(hex_grid_path(cfg), cons_areas_hex_path(cfg)))),
    resolution = cfg$hex_resolution,
    crs        = cfg$crs$working,
    code       = lapply(fns, function(f) deparse(get(f, envir = globalenv())))
  ))
}

### Returns NULL if the taxon's outputs are reusable, else a short reason string.
range_stale_reason <- function(taxon, cfg = config, taxa_cfg = taxa_config) {
  p <- range_manifest_path(taxon, cfg)
  if (!file.exists(p)) return("no manifest")
  m <- tryCatch(jsonlite::read_json(p, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(m)) return("unreadable manifest")
  outs <- .range_output_paths(taxon, cfg)
  if (!all(file.exists(outs))) return("output file missing")
  rec <- unlist(m$outputs$md5)
  if (!identical(unname(rec[names(outs)]), unname(tools::md5sum(outs)))) {
    return("output changed on disk")
  }
  if (!identical(m$settings, range_settings_hash(taxon, cfg, taxa_cfg))) return("inputs changed")
  NULL
}

.write_range_manifest <- function(taxon, cfg, taxa_cfg) {
  outs <- .range_output_paths(taxon, cfg)
  rc <- arrow::open_dataset(outs[["range_cells"]])
  mpath <- range_manifest_path(taxon, cfg)
  ensure_parent_dir(mpath)
  jsonlite::write_json(
    list(taxon = taxon, range_source = taxa_cfg[[taxon]]$range_source,
         settings = range_settings_hash(taxon, cfg, taxa_cfg),
         n_cells = nrow(rc), n_species = ncol(rc) - 1L,
         outputs = list(paths = as.list(sub(paste0("^", here::here(), "/"), "", outs)),
                        md5 = as.list(setNames(unname(tools::md5sum(outs)), names(outs)))),
         built_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S")),
    mpath, auto_unbox = TRUE, pretty = TRUE
  )
}

run_build_ranges <- function(taxa = config$taxa, cfg = config, taxa_cfg = taxa_config,
                             force = FALSE,
                             progress = c("terse", "verbose", "quiet")) {
  progress <- match.arg(progress)

  ### loop over taxa
  for (taxon in taxa) {
    ### access the range source from the taxa_config:
    src <- taxa_cfg[[taxon]]$range_source

    ### check that source is valid:
    if (!is_polygon_source(src) && !identical(src, "raster")) {
      warning(glue::glue("[{taxon}] unknown range_source '{src}'; skipped."))
      next
    }

    ### if force == FALSE, check if cache is stale; if not, next
    if (!isTRUE(force)) {
      reason <- range_stale_reason(taxon, cfg, taxa_cfg)
      if (is.null(reason)) {
        .progress_msg(progress, "terse",
                      glue::glue("[range] {taxon}: CACHED (inputs unchanged; force = TRUE to rebuild)"))
        next
      }
      .progress_msg(progress, "terse", glue::glue("[range] {taxon} (source: {src}; {reason})"))
    } else {
      .progress_msg(progress, "terse", glue::glue("[range] {taxon} (source: {src}; forced)"))
    }
    ### Drop any old manifest first so a failed rebuild can't leave a stale "valid" one.
    unlink(range_manifest_path(taxon, cfg))

    ### build the range from the source as appropriate
    if (is_polygon_source(src)) {
      build_range_from_polygons(taxon, cfg, taxa_cfg, progress = progress)
    } else {
      build_range_from_raster(taxon, cfg, taxa_cfg, progress = progress)
    }
    .write_range_manifest(taxon, cfg, taxa_cfg)
  }
  invisible(taxa)
}
