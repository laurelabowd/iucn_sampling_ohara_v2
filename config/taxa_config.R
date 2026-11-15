### =============================================================================
### config/taxa_config.R  —  per-taxon settings
### =============================================================================
###
### One entry per taxon. Downstream code never hardcodes "bird"/"reptile"/"mammal"
### or a per-taxon path — it loops over `names(taxa_config)` and reads from here.
###
### Requires config/config.R to be sourced first (uses `config`).
###
### Each taxon defines:
###   species_source        where the raw species list comes from (CSV or Google
###                         Sheet) and which column holds the scientific name.
###   range_source          "polygon" -> ranges come from per-species range
###                         vector source (e.g., IUCN rangemap vector data).
###                         "raster"  -> ranges come from a raster source (e.g.,
###                         eBird data or AquaX).  Note, "shapefile" or "vector" 
###                         also work for vector/polygon data
###   polygon_dir           (polygon taxa) folder of per-species range polygons.
###   polygon_name_field    (polygon taxa) attribute holding the species name in
###                         each polygon (the source's own taxonomy; Stage 0
###                         resolves it to the GBIF key so it can be joined).
###   obis_dataset_path     (marine taxa) path to OBIS occurrence dataset (parquet).
###   gbif_dataset_path     raw GBIF/OBIS occurrence dataset (parquet) for this taxon,
###                         parsed into per-species-year files in 03_process.

library(here)

### Columns kept from raw GBIF when parsing. Defined ONCE and shared by every
### taxon so a per-taxon list can't silently drop a column the cleaning step needs
### (the original mammal config dropped basisOfRecord/individualCount, which
### clean_occurrence then filtered on — a real bug this shared list prevents).
gbif_columns_keep <- c(
  "gbifID", "datasetKey", "occurrenceID",
  "family", "genus", "species", "scientificName",
  "occurrenceStatus", "publishingOrgKey",
  "decimalLatitude", "decimalLongitude",
  "coordinateUncertaintyInMeters", "coordinatePrecision",
  "year", "month", "day", "eventDate",
  "basisOfRecord", "individualCount"
)

### Helper: default species source is the shared sample list, filtered to `taxon`
### (the sample CSV has a `taxon` column). Point csv_path at your own file, or set
### sheet_url/sheet_name to read a Google Sheet tab instead.
.species_source <- function(taxon, sci_name_col = "scientific_name", csv_path = NULL) {
  list(
    csv_path     = if (is.null(csv_path)) config$species_list_default else csv_path,
    taxon_filter = taxon,          ### filter csv rows to this taxon (NULL to skip)
    sheet_url    = NULL,
    sheet_name   = NULL,
    sci_name_col = sci_name_col
  )
}

taxa_config <- list(

  # bird = list(
  #   taxon                = "bird",
  #   species_source       = .species_source("bird"),
  #   range_source         = "raster",
  #   ebird_richness_raster = file.path(config$paths$data_processed,
  #                                     "species_populations", "bird",
  #                                     "ebird_range_richness.tif"),
  #   gbif_dataset_path    = file.path(config$paths$data_raw, "gbif", "bird",
  #                                    "occurrence.parquet"),
  #   gbif_columns_keep    = gbif_columns_keep
  # ),

  ### Marine mammals: species list from IUCN spatial data extraction
  ### see R/_pre_analysis scripts, ranges from IUCN Red List polygons
  mammals_marine = list(
    taxon                = "mammals_marine",
    species_source       = .species_source(
      "mammals_marine",
      csv_path = file.path(config$paths$data_raw, "iucn_species_list.csv")
    ),
    range_source         = "polygon",
    polygon_dir          = file.path(config$paths$data_raw, "iucn", "mammals_marine"),
    polygon_name_field   = "sci_name",    ### IUCN field holding species name
    obis_dataset_path    = file.path(config$paths$data_raw, "obis", "mammals_marine",
                                     "occurrence.parquet"),
    gbif_dataset_path    = file.path(config$paths$data_raw, "gbif", "mammals_marine",
                                     "occurrence.parquet"),
    gbif_columns_keep    = gbif_columns_keep
  ),

  ### Abalones: species list from IUCN spatial data extraction
  abalones = list(
    taxon                = "abalones",
    species_source       = .species_source(
      "abalones",
      csv_path = file.path(config$paths$data_raw, "iucn_species_list.csv")
    ),
    range_source         = "polygon",
    polygon_dir          = file.path(config$paths$data_raw, "iucn", "abalones"),
    polygon_name_field   = "sci_name",    ### IUCN field holding species name
    obis_dataset_path    = file.path(config$paths$data_raw, "obis", "abalones",
                                     "occurrence.parquet"),
    gbif_dataset_path    = file.path(config$paths$data_raw, "gbif", "abalones",
                                     "occurrence.parquet"),
    gbif_columns_keep    = gbif_columns_keep
  )
)
