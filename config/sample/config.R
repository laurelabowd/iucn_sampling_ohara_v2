### =============================================================================
### config/sample/config.R  —  global configuration for the OFFLINE SAMPLE run
### =============================================================================
###
### Sourced by apply_sample_config() in R/run_pipeline.R, which overwrites the
### in-memory `config` from config/config.R. All inputs point into data/sample/
### and field names/values match what data/sample/make_sample_data.R wrote.
### No API calls: the species crosswalk is committed.

library(here)

config <- list(

  ### --- Study region ----------------------------------------------------------
  region = list(
    country_code    = NULL,
    admin_level     = NULL,
    admin_name      = 'sample_region',
    label           = 'sample_region',
    boundary_source = here('data', 'sample', 'spatial_layers', 'region.gpkg'),
    bbox            = c(-120, -119, 35, 36)
  ),

  ### --- Conservation-areas layer (single-layer gpkg, sample schema) -----------
  conservation_areas = list(
    path               = here('data', 'sample', 'spatial_layers', 'conservation_areas.gpkg'),
    layer              = NULL,
    gap_field          = 'reGAP',
    gap_values         = c(1, 2),
    realm_field        = 'TerrMar',
    realm_value        = 'Terrestrial',
    unit_field         = 'UNIT_NAME',
    agency_field       = 'MNG_AGNCY',
    agency_level_field = 'MNG_AG_LEV'
  ),

  ### --- Years -----------------------------------------------------------------
  years         = c(2014, 2024),
  focal_year    = 2024,
  compare_years = c(2014, 2024),

  ### --- Coordinate reference systems ------------------------------------------
  crs = list(
    working    = "EPSG:4326",   # H3 / geographic
    equal_area = "EPSG:3310",   # California Albers — swap for your region
    display    = "EPSG:3857"    # web mercator, plotting only
  ),

  ### --- Spatial resolution ----------------------------------------------------
  hex_resolution = 8,

  ### --- API credentials (unused offline; kept so downstream code finds the keys)
  credentials = list(
    iucn_redlist_key = Sys.getenv('IUCN_REDLIST_KEY'),
    gbif_user        = Sys.getenv('GBIF_USER'),
    gbif_pwd         = Sys.getenv('GBIF_PWD'),
    gbif_email       = Sys.getenv('GBIF_EMAIL')
  ),

  ### --- Occurrence sources: offline sample uses the committed GBIF parquet only
  gbif_filters = list(
    enabled         = TRUE,
    has_geo_issue   = FALSE,
    has_coordinate  = TRUE,
    basis_of_record = c('HUMAN_OBSERVATION', 'OBSERVATION',
                        'MACHINE_OBSERVATION', 'PRESERVED_SPECIMEN'),
    download_format = 'SIMPLE_PARQUET'
  ),
  obis_filters = list(
    enabled         = FALSE,
    basis_of_record = c('HUMAN_OBSERVATION', 'OBSERVATION',
                        'MACHINE_OBSERVATION', 'PRESERVED_SPECIMEN')
  ),
  obis_use_aws      = FALSE,
  obis_max_age_days = 180,
  force             = TRUE,
  gbif_max_age_days = 180,
  ### 'download' = GBIF download API, so the API path is exercised on sample runs
  gbif_source       = 'download',
  ### GBIF backbone checklist key to use GBIF Backbone Taxonomy for name matching
  gbif_checklist_key = 'd7dddbf4-2cf0-4f39-9b2a-bb099caae36c',

  ### --- Paths (sample-scoped so nothing collides with real outputs) ----------
  paths = list(
    data_raw       = here('data', 'sample'),
    data_processed = here('data', 'sample', 'processed'),
    output         = here('data', 'sample', 'output'),
    crosswalk      = here('data', 'sample', 'species_crosswalk.csv'),
    hex_grid       = here('data', 'sample', 'processed', 'spatial_layers', 'hex_grid.gpkg'),
    cons_areas_hex = here('data', 'sample', 'processed', 'spatial_layers', 'cons_areas_hex.gpkg')
  ),

  ### --- Occurrence-cleaning thresholds (same as full run) --------------------
  cleaning = list(
    max_coord_uncertainty_m = 750,
    min_coord_precision     = 0.01,
    max_individual_count    = 99,
    exclude_basis_of_record = c('FOSSIL_SPECIMEN', 'LIVING_SPECIMEN')
  ),

  ### --- Species list ----------------------------------------------------------
  species_list_default = here('data', 'sample', 'species_list.csv'),

  ### --- Taxa (must match names in config/sample/taxa_config.R) ---------------
  taxa = c('reptile', 'mammal')
)

### Output-path helpers used by the stages; mirror config/config.R exactly.
taxon_output_dir <- function(taxon, cfg = config) {
  file.path(cfg$paths$output, 'indicators', 'species_populations', taxon)
}

taxon_processed_dir <- function(taxon, cfg = config) {
  file.path(cfg$paths$data_processed, 'species_populations', taxon)
}

taxon_year_occ_dir <- function(taxon, year, cfg = config) {
  file.path(taxon_processed_dir(taxon, cfg), 'occurrences', as.character(year))
}
