### =============================================================================
### config/sample/taxa_config.R  —  per-taxon settings for the OFFLINE SAMPLE run
### =============================================================================
###
### Two polygon taxa (reptile, mammal) built from data/sample/. Species-to-range
### links come from data/sample/species_crosswalk.csv (range_polygon_path), and
### each range shapefile stores the species name in `SName`.
###
### Requires config/sample/config.R sourced first (uses `config`).

library(here)

### Same column set as config/taxa_config.R so the parser sees the fields it expects.
gbif_columns_keep <- c(
  "gbifID", "datasetKey", "occurrenceID",
  "family", "genus", "species", "scientificName",
  "occurrenceStatus", "publishingOrgKey",
  "decimalLatitude", "decimalLongitude",
  "coordinateUncertaintyInMeters", "coordinatePrecision",
  "year", "month", "day", "eventDate",
  "basisOfRecord", "individualCount"
)

taxa_config <- list(

  reptile = list(
    taxon              = "reptile",
    species_source     = list(
      csv_path     = config$species_list_default,
      taxon_filter = "reptile",
      sheet_url    = NULL,
      sheet_name   = NULL,
      sci_name_col = "scientific_name"
    ),
    range_source       = "polygon",
    polygon_dir        = here("data", "sample", "ranges", "reptile"),
    polygon_name_field = "SName",
    gbif_dataset_path  = here("data", "sample", "gbif", "reptile", "occurrence.parquet"),
    gbif_columns_keep  = gbif_columns_keep
  ),

  mammal = list(
    taxon              = "mammal",
    species_source     = list(
      csv_path     = config$species_list_default,
      taxon_filter = "mammal",
      sheet_url    = NULL,
      sheet_name   = NULL,
      sci_name_col = "scientific_name"
    ),
    range_source       = "polygon",
    polygon_dir        = here("data", "sample", "ranges", "mammal"),
    polygon_name_field = "SName",
    gbif_dataset_path  = here("data", "sample", "gbif", "mammal", "occurrence.parquet"),
    gbif_columns_keep  = gbif_columns_keep
  )
)
