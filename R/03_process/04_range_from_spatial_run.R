### =============================================================================
### R/03_process/04_range_from_spatial_run.R  —  range-building driver
### =============================================================================
### Builds per-species hex range cells for every taxon, dispatching by
### range_source: "polygon" -> 04a_range_from_polygons_fn.R (aliases: "shapefile", "vector")
###               "raster"  -> 04b_range_from_raster_fn.R (eBird birds).
### `progress`: "terse" (default), "verbose" (per layer/species detail incl. area),
### or "quiet". Both write the same range_grid_cells.parquet shape, so downstream is uniform.
###
### Requires config + taxa_config + sampling_config + R/00_utils.R.

source(here::here("R", "03_process", "04a_range_from_polygons_fn.R"))
source(here::here("R", "03_process", "04b_range_from_raster_fn.R"))

run_build_ranges <- function(taxa = config$taxa, cfg = config, taxa_cfg = taxa_config,
                             progress = c("terse", "verbose", "quiet")) {
  progress <- match.arg(progress)
  for (taxon in taxa) {
    src <- taxa_cfg[[taxon]]$range_source
    .progress_msg(progress, "terse", glue::glue("[range] {taxon} (source: {src})"))
    if (is_polygon_source(src)) {
      build_range_from_polygons(taxon, cfg, taxa_cfg, progress = progress)
    } else if (identical(src, "raster")) {
      build_range_from_raster(taxon, cfg, taxa_cfg, progress = progress)
    } else {
      warning(glue::glue("[{taxon}] unknown range_source '{src}'; skipped."))
    }
  }
  invisible(taxa)
}
