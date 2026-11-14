### =============================================================================
### R/03_process/04_range_from_shapefiles_run.R  —  range-building driver
### =============================================================================
### Builds per-species hex range cells for every taxon, dispatching by
### range_source: "shapefile" -> 04a_range_from_shapefiles_fn.R,
###               "raster"    -> 04b_range_from_raster_fn.R (eBird birds).
### Both write the same range_grid_cells.parquet shape, so downstream is uniform.
###
### Requires config + taxa_config + sampling_config + R/00_utils.R.

source(here::here("R", "03_process", "04a_range_from_shapefiles_fn.R"))
source(here::here("R", "03_process", "04b_range_from_raster_fn.R"))

run_build_ranges <- function(taxa = config$taxa, cfg = config, taxa_cfg = taxa_config) {
  for (taxon in taxa) {
    src <- taxa_cfg[[taxon]]$range_source
    message(glue::glue("[range] {taxon} (source: {src}) ..."))
    if (identical(src, "shapefile")) {
      build_range_from_shapefiles(taxon, cfg, taxa_cfg)
    } else if (identical(src, "raster")) {
      build_range_from_raster(taxon, cfg, taxa_cfg)
    } else {
      warning(glue::glue("[{taxon}] unknown range_source '{src}'; skipped."))
    }
  }
  invisible(taxa)
}
