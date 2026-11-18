### =============================================================================
### config/sample/sampling_config.R  —  sampling-optimization config, OFFLINE SAMPLE
### =============================================================================
###
### Same tuning knobs as config/sampling_config.R. Requires config/sample/config.R
### and config/sample/taxa_config.R sourced first.

sampling_params <- list(
  well_sampled_threshold = 50,
  target_sampling_frac   = 0.10,
  cons_balance_weight    = 0.20,
  batch_fraction         = 0.35,
  deficit_percentile     = 0.15,
  hex_resolution         = config$hex_resolution
)

### Same builder as config/sampling_config.R; defaults point at the sample config.
make_sampling_config <- function(year, suffix = "", cfg = config,
                                 taxa_cfg = taxa_config) {
  year <- as.character(year)
  tag  <- if (nzchar(suffix)) paste0("_", suffix) else ""

  build_one <- function(taxon) {
    tcfg     <- taxa_cfg[[taxon]]
    proc_dir <- taxon_processed_dir(taxon, cfg)
    out_dir  <- taxon_output_dir(taxon, cfg)

    list(
      taxon        = taxon,
      range_source = tcfg$range_source,
      occ_dir      = taxon_year_occ_dir(taxon, year, cfg),
      range_cells_path = file.path(proc_dir, "range_grid_cells.parquet"),
      hex_richness_path = file.path(out_dir,
                                    paste0("hex_species_richness_", year, ".gpkg")),
      cell_list_path    = file.path(out_dir, "spp_cells_in_region.csv"),
      output_dir        = file.path(out_dir, "sampling_optimized"),
      undersamp_cells_out = file.path(proc_dir,
                                      paste0("undersamp_cells_", year, tag, ".parquet")),
      undersamp_spp_out   = file.path(proc_dir,
                                      paste0("undersamp_spp_", year, tag, ".csv"))
    )
  }

  stats::setNames(lapply(cfg$taxa, build_one), cfg$taxa)
}
