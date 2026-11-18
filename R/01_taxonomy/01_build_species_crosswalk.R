### =============================================================================
### R/01_taxonomy/01_build_species_crosswalk.R  —  STAGE 0 driver
### =============================================================================
###
### Turns each taxon's raw input species list into ONE authoritative crosswalk
### that every downstream stage joins on. Runs FIRST in the pipeline.
###
### What it adds beyond resolve_species_names():
###   - range_polygon_name / range_polygon_path / range_polygon_layer for polygon-range taxa.
###     Range polygons carry the vendor's own species name (the polygon_name_
###     field, e.g. "SName"), which may not match the GBIF canonical name. We
###     resolve those polygon names through the GBIF backbone and join them to
###     the crosswalk on the shared gbif_usage_key -- the same key-bridge the
###     resolver already uses for eBird. This is the join key that lets range
###     cells (keyed by polygon name) line up with occurrence/richness data
###     (keyed by GBIF canonical name).
###
###   Polygon sources may be multi-species (e.g. IUCN geopackages holding many
###   species as rows of one layer) or single-species (one file per species), so
###   the link records the file path, the layer within it, and the species name.
###   Rows of a layer are later selected by matching the name field.
###
### Output: data/processed/species_crosswalk.csv (one row per species, all taxa).
###
### Requires config/config.R + config/taxa_config.R sourced first.
###
### Usage:
###   source("config/config.R"); source("config/taxa_config.R")
###   source("R/00_utils.R")
###   source("R/01_taxonomy/01_build_species_crosswalk.R")
###   crosswalk <- build_species_crosswalk()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(glue)
  library(janitor)
  library(sf)
  library(here)
})

source(here::here("R", "01_taxonomy", "01a_resolve_species_names.R"))

### Load one taxon's raw species list from its configured CSV or Google Sheet.
### When reading a shared CSV that carries a `taxon` column, filter to this taxon.
load_species_list <- function(taxon, tcfg) {
  src <- tcfg$species_source
  has_sheet <- !is.null(src$sheet_url) && !is.na(src$sheet_url) && nzchar(src$sheet_url)

  if (has_sheet) {
    googlesheets4::gs4_deauth()
    df <- googlesheets4::read_sheet(ss = src$sheet_url, sheet = src$sheet_name)
  } else if (!is.null(src$csv_path) && file.exists(src$csv_path)) {
    df <- readr::read_csv(src$csv_path, show_col_types = FALSE)
  } else {
    stop(glue("No species source found for taxon '{taxon}' (csv_path missing)."))
  }

  ### Optional filter when a shared list holds multiple taxa.
  if (!is.null(src$taxon_filter)) {
    df_clean <- janitor::clean_names(df)
    if ("taxon" %in% names(df_clean)) {
      df <- df[tolower(df_clean$taxon) == tolower(src$taxon_filter), , drop = FALSE]
    }
  }
  df
}

### Count features per species name in one layer of a vector file, attributes
### only. Tries an OGR SQL query first (avoids loading geometry for big layers),
### falling back to a full read. Returns tibble(name, n_features) or NULL.
.layer_polygon_names <- function(path, layer, name_field) {
  q <- glue('SELECT "{name_field}" AS nm, COUNT(*) AS n FROM "{layer}" GROUP BY "{name_field}"')
  counts <- tryCatch({
    res <- suppressWarnings(sf::st_read(path, query = q, quiet = TRUE))
    tibble::tibble(name = as.character(res$nm), n_features = as.integer(res$n))
  }, error = function(e) NULL)

  if (is.null(counts)) {
    atts <- sf::st_read(path, layer = layer, quiet = TRUE) |> sf::st_drop_geometry()
    if (!name_field %in% names(atts)) {
      warning(glue("'{basename(path)}' layer '{layer}' has no field '{name_field}'; skipped"))
      return(NULL)
    }
    counts <- dplyr::count(tibble::tibble(name = as.character(atts[[name_field]])),
                           name, name = "n_features")
  }
  dplyr::filter(counts, !is.na(name), name != "")
}

### Read the species-name attribute out of every vector file in a directory.
### Handles multi-species files (many species per layer, e.g. IUCN geopackages)
### and single-species files alike: every layer of every .gpkg, plus each .shp
### (whose layer is its file stem). Returns one row per species name per layer:
### range_polygon_name, range_polygon_path, range_polygon_layer, range_polygon_n_features.
harvest_polygon_names <- function(polygon_dir, name_field) {
  empty <- tibble::tibble(range_polygon_name = character(),
                          range_polygon_path = character(),
                          range_polygon_layer = character(),
                          range_polygon_n_features = integer())
  if (is.null(polygon_dir) || !dir.exists(polygon_dir)) return(empty)

  files <- list.files(polygon_dir, pattern = "\\.(shp|gpkg)$",
                      full.names = TRUE, recursive = TRUE)
  purrr::map_dfr(files, function(f) {
    tryCatch({
      layers <- if (grepl("\\.gpkg$", f)) sf::st_layers(f)$name
                else tools::file_path_sans_ext(basename(f))
      purrr::map_dfr(layers, function(lyr) {
        cnt <- .layer_polygon_names(f, lyr, name_field)
        if (is.null(cnt) || nrow(cnt) == 0) return(NULL)
        tibble::tibble(range_polygon_name = cnt$name, range_polygon_path = f,
                       range_polygon_layer = lyr, range_polygon_n_features = cnt$n_features)
      })
    }, error = function(e) {
      warning(glue("could not read '{basename(f)}': {conditionMessage(e)}"))
      NULL
    })
  }) |>
    dplyr::bind_rows(empty)
}

### Resolve polygon names to GBIF keys and attach them to the crosswalk on
### gbif_usage_key. Rows in the crosswalk with no matching polygon get NA.
attach_polygon_link <- function(crosswalk, polygon_names) {
  if (nrow(polygon_names) == 0) {
    return(dplyr::mutate(crosswalk,
                         range_polygon_name = NA_character_,
                         range_polygon_path = NA_character_,
                         range_polygon_layer = NA_character_,
                         range_polygon_n_features = NA_integer_))
  }

  keys <- rgbif::name_backbone_checklist(polygon_names$range_polygon_name) |>
    dplyr::transmute(
      range_polygon_name = verbatim_name,
      gbif_usage_key = as.numeric(speciesKey)
    )

  ### One crosswalk row per species. If a species has polygons in several
  ### layers/files (or under several vendor names), the sources are collapsed into
  ### parallel "|"-delimited name/path/layer strings that the range step expands
  ### and combines. n_features totals polygon features across all sources.
  link <- polygon_names |>
    dplyr::left_join(keys, by = "range_polygon_name") |>
    dplyr::filter(!is.na(gbif_usage_key)) |>
    dplyr::summarise(
      range_polygon_name = paste(range_polygon_name, collapse = "|"),
      range_polygon_path = paste(range_polygon_path, collapse = "|"),
      range_polygon_layer = paste(range_polygon_layer, collapse = "|"),
      range_polygon_n_features = sum(range_polygon_n_features),
      .by = gbif_usage_key
    )

  multi <- dplyr::filter(link, range_polygon_n_features > 1)
  if (nrow(multi) > 0) {
    warning(glue("{nrow(multi)} species have >1 polygon feature (e.g. differing ",
                 "presence/seasonality); the range step will combine them: ",
                 "{paste(utils::head(multi$range_polygon_name, 5), collapse = ', ')}",
                 "{if (nrow(multi) > 5) ', ...' else ''}"), call. = FALSE)
  }

  crosswalk |>
    dplyr::left_join(link, by = "gbif_usage_key")
}

### Build one taxon's crosswalk (names resolved + polygon linkage where used).
build_taxon_crosswalk <- function(taxon, tcfg, cfg = config, force = config$force) {
  raw <- load_species_list(taxon, tcfg)

  cw <- resolve_species_names(
    spp_df = raw,
    sci_name_col = tcfg$species_source$sci_name_col,
    taxon_group = taxon,          ### stamp group so reptiles aren't left NA
    rebuild_ebird_cache = isTRUE(force)
  )

  if (is_polygon_source(tcfg$range_source)) {
    polygon_names <- harvest_polygon_names(tcfg$polygon_dir, tcfg$polygon_name_field)
    cw <- attach_polygon_link(cw, polygon_names)

    ### Flag species whose range polygon we could not match -> review.
    cw <- cw |>
      dplyr::mutate(
        needs_review = needs_review | is.na(range_polygon_name),
        review_reason = dplyr::if_else(
          is.na(range_polygon_name) & !override,
          stringr::str_trim(paste(review_reason, "no range polygon matched", sep = "; ")),
          review_reason
        )
      )
  } else {
    cw <- dplyr::mutate(cw,
                        range_polygon_name = NA_character_,
                        range_polygon_path = NA_character_,
                        range_polygon_layer = NA_character_,
                        range_polygon_n_features = NA_integer_)
  }

  dplyr::mutate(cw, taxon = taxon, .before = 1)
}

### Build the master crosswalk across all configured taxa and write it once.
### `force = TRUE` also refreshes the cached eBird-to-GBIF key map.
build_species_crosswalk <- function(taxa = config$taxa,
                                    taxa_cfg = taxa_config,
                                    cfg = config,
                                    write = TRUE,
                                    force = config$force) {
  ### test to warn if config and taxa_config don't match for all taxa!
  all_missed <- c(setdiff(taxa, names(taxa_config)), setdiff(names(taxa_config), taxa))
  if(length(all_missed) != 0) {
    warning('Unmatched taxa!',
            '\n  Taxa in config:      ', paste0(taxa, collapse = ', '),
            '\n  Taxa in taxa_config: ', paste0(names(taxa_cfg), collapse = ', '))
  }
  
  crosswalk <- purrr::map_dfr(taxa, function(taxon) {
    tcfg <- taxa_cfg[[taxon]]
    if (is.null(tcfg)) stop(glue("Taxon '{taxon}' not found in taxa_config."))
    build_taxon_crosswalk(taxon, tcfg, cfg, force = force)
  })

  if (isTRUE(write)) {
    out <- cfg$paths$crosswalk
    ensure_parent_dir(out)
    readr::write_csv(crosswalk, out)
    message(glue(
      "[crosswalk] {nrow(crosswalk)} species across {length(taxa)} taxa -> {out} ",
      "({sum(crosswalk$needs_review)} flagged for review)"
    ))
  }

  crosswalk
}
