### =============================================================================
### R/02_acquire/02_download_obis.R  —  download OBIS occurrence records
### =============================================================================
###
### Marine counterpart of 01_download_gbif.R. For large taxon downloads,
### queries the OBIS GeoParquet dataset hosted on AWS S3 (via DuckDB) instead
### of paging through the API. For small downloads or when S3 access is unavailable,
### falls back to the OBIS API.
###
### Large downloads (>100k records estimated) use:
###   - AWS S3 GeoParquet files at s3://obis-open-data/occurrence/*.parquet
###   - DuckDB for spatial + attribute filtering (faster than API paging)
###   - Direct parquet write to taxa_config$obis_dataset_path
###
### Small downloads use the OBIS API (https://api.obis.org/v3/occurrence) with
### species-name matching, paged through 10k-record batches.
###
### Filters apply uniformly: year range (config$years), region bounding box,
### and basis_of_record (config$obis_filters). Columns are renamed to GBIF-style
### names so parse/clean stages treat both sources alike.
###
### No credentials needed. Skipped entirely when config$obis_filters$enabled is FALSE.
###
### Requires config + taxa_config + R/00_utils.R + R/02_acquire/01_download_gbif.R
### (for load_crosswalk() and region_bbox()) sourced first.
###
### Usage:
###   res <- run_obis_downloads()                 # all configured taxa
###   res <- run_obis_downloads(taxa = "mammals_marine", force = TRUE)
###   res <- run_obis_downloads(use_aws = FALSE)  # API only, no S3

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(glue)
  library(httr2)
  library(arrow)
  library(here)
  library(DBI)
})

obis_api_url <- "https://api.obis.org/v3/occurrence"
obis_s3_prefix <- "s3://obis-open-data/occurrence"

### Fields requested from OBIS (these are also the columns written to parquet).
obis_fields <- c(
  "id", "datasetID", "aphiaID", "scientificName", "species",
  "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
  "date_year", "month", "day", "eventDate",
  "basisOfRecord", "individualCount", "occurrenceStatus"
)

### Accepted binomials to query. Flagged (needs_review) rows are included by
### default; set include_flagged = FALSE to drop them.
get_obis_names <- function(crosswalk, include_flagged = TRUE) {
  nms <- dplyr::mutate(crosswalk, query_name = dplyr::coalesce(name_gbif, name_accepted))
  nms <- dplyr::filter(nms, !is.na(query_name))
  if (!isTRUE(include_flagged)) nms <- dplyr::filter(nms, !needs_review)
  sort(unique(nms$query_name))
}

### WKT polygon for the region, or NULL when no spatial restriction is needed.
### Reuses region_bbox() from 01_download_gbif.R.
obis_region_wkt <- function(cfg = config) {
  bb <- region_bbox(cfg)
  if (is.null(bb)) return(NULL)
  as.character(glue::glue(
    "POLYGON(({bb[['xmin']]} {bb[['ymin']]}, {bb[['xmax']]} {bb[['ymin']]}, ",
    "{bb[['xmax']]} {bb[['ymax']]}, {bb[['xmin']]} {bb[['ymax']]}, ",
    "{bb[['xmin']]} {bb[['ymin']]}))"
  ))
}

### --- AWS S3 GeoParquet path (new, faster path) -------------------------------------------
### Query OBIS GeoParquet on AWS S3 via DuckDB. Much faster for large datasets than API paging.
### Returns a data frame with OBIS fields (or NULL if DuckDB unavailable or query fails).
fetch_obis_aws <- function(species_names, year_range, bbox = NULL, cfg = config) {
  # DuckDB support is optional; return NULL if not available
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    db <- duckdb::dbConnect(duckdb::duckdb(), read_only = TRUE)
    on.exit(duckdb::dbDisconnect(db), add = TRUE)
    
    # Enable S3 access (auto-install and auto-load httpfs extension)
    DBI::dbExecute(db, "SET autoinstall_known_extensions=1")
    DBI::dbExecute(db, "SET autoload_known_extensions=1")
    
    ### The GeoParquet nests QC'd fields in the `interpreted` struct (see
    ### update_docs/08_AWS_OBIS_SCHEMA_DEBUG.md). Alias to the names in
    ### obis_fields so standardize_obis() treats S3 and API results alike.
    bbox_filter <- if (!is.null(bbox)) {
      glue::glue(
        "AND interpreted.decimalLongitude >= {bbox[['xmin']]} ",
        "AND interpreted.decimalLongitude <= {bbox[['xmax']]} ",
        "AND interpreted.decimalLatitude >= {bbox[['ymin']]} ",
        "AND interpreted.decimalLatitude <= {bbox[['ymax']]}"
      )
    } else ""

    ### dbQuoteString escapes names containing apostrophes (shQuote does not).
    name_list <- paste(DBI::dbQuoteString(db, species_names), collapse = ", ")

    sql <- glue::glue("
      SELECT
        CAST(_id AS VARCHAR) AS id,
        CAST(dataset_id AS VARCHAR) AS datasetID,
        interpreted.aphiaid AS aphiaID,
        interpreted.scientificName AS scientificName,
        interpreted.species AS species,
        interpreted.decimalLatitude AS decimalLatitude,
        interpreted.decimalLongitude AS decimalLongitude,
        interpreted.coordinateUncertaintyInMeters AS coordinateUncertaintyInMeters,
        interpreted.date_year AS date_year,
        interpreted.month AS month,
        interpreted.day AS day,
        interpreted.eventDate AS eventDate,
        interpreted.basisOfRecord AS basisOfRecord,
        interpreted.individualCount AS individualCount,
        interpreted.occurrenceStatus AS occurrenceStatus
      FROM read_parquet('{obis_s3_prefix}/*.parquet')
      WHERE
        interpreted.species IN ({name_list})
        AND interpreted.date_year >= {as.integer(year_range[1])}
        AND interpreted.date_year <= {as.integer(year_range[2])}
        {bbox_filter}
    ")
    
    # Use DBI interface: dbGetQuery or dbExecute + dbFetch
    recs <- DBI::dbGetQuery(db, sql)
    dplyr::as_tibble(recs)
  },
  error = function(e) {
    warning(glue("AWS S3 query failed: {conditionMessage(e)}; falling back to API."))
    NULL
  })
}

### --- Download manifest -------------------------------------------------------
### Same idea as the GBIF manifest: records what was downloaded and from which
### request, so unchanged requests are not re-pulled.
obis_manifest_path <- function(out_dir) file.path(out_dir, "obis_download_manifest.json")

### Hash of everything that determines the download's content.
obis_request_hash <- function(obis_names, cfg = config) {
  rlang::hash(list(
    names   = obis_names,
    filters = cfg$obis_filters$basis_of_record,
    years   = range(cfg$years),
    region  = region_bbox(cfg),
    fields  = obis_fields
  ))
}

read_obis_manifest <- function(out_dir) {
  p <- obis_manifest_path(out_dir)
  if (!file.exists(p)) return(NULL)
  tryCatch(jsonlite::read_json(p, simplifyVector = TRUE), error = function(e) NULL)
}

### Returns NULL if the existing download is reusable, else a short reason string.
obis_download_stale_reason <- function(manifest, request_hash, out_file,
                                       max_age_days = Inf) {
  if (is.null(manifest)) return("no manifest")
  if (!identical(manifest$request_hash, request_hash)) return("request changed")
  if (!file.exists(out_file)) return("downloaded file missing")
  age <- as.numeric(difftime(Sys.time(), as.POSIXct(manifest$downloaded_utc, tz = "UTC"),
                             units = "days"))
  if (is.na(age) || age > max_age_days) {
    return(glue("older than max_age_days ({round(age, 1)} d)"))
  }
  NULL
}

### Rename to GBIF-style columns and apply the shared basis-of-record filter.
### OBIS uses CamelCase ("HumanObservation"); GBIF uses "HUMAN_OBSERVATION".
standardize_obis <- function(df, cfg = config) {
  if (nrow(df) == 0) return(df)
  ### Ensure every expected column exists even if no record populated it.
  for (col in setdiff(obis_fields, names(df))) df[[col]] <- NA
  df <- df |>
    dplyr::select(dplyr::all_of(obis_fields)) |>
    dplyr::rename(year = date_year, obisID = id) |>
    dplyr::mutate(
      basisOfRecord = toupper(gsub("([a-z])([A-Z])", "\\1_\\2", basisOfRecord)),
      year = as.integer(year)
    )
  keep <- cfg$obis_filters$basis_of_record
  if (!is.null(keep)) df <- dplyr::filter(df, basisOfRecord %in% keep)
  dplyr::distinct(df, obisID, .keep_all = TRUE)
}

### Page through the API for one batch of names. OBIS pages by an `after` cursor
### (the last record id), so we loop until an empty page comes back.
fetch_obis_batch <- function(names_batch, year_range, wkt = NULL, page_size = 10000) {
  query <- list(
    scientificname = paste(names_batch, collapse = ","),
    startdate = sprintf("%d-01-01", as.integer(year_range[[1]])),
    enddate   = sprintf("%d-12-31", as.integer(year_range[[2]])),
    size      = page_size,
    fields    = paste(obis_fields, collapse = ",")
  )
  if (!is.null(wkt)) query$geometry <- wkt

  pages <- list()
  after <- NULL
  repeat {
    q <- query
    if (!is.null(after)) q$after <- after
    body <- httr2::request(obis_api_url) |>
      httr2::req_url_query(!!!q) |>
      httr2::req_retry(max_tries = 4) |>
      httr2::req_perform() |>
      httr2::resp_body_json()
    res <- body$results
    if (length(res) == 0) break
    pages[[length(pages) + 1]] <- dplyr::bind_rows(lapply(res, tibble::as_tibble_row))
    after <- res[[length(res)]]$id
    if (length(res) < page_size) break
  }
  dplyr::bind_rows(pages)
}

### Download one taxon; returns a one-row result tibble (a run-log entry).
### Skips when the manifest matches the current request and is fresh (see
### `max_age_days`); `force = TRUE` always re-downloads.
### use_aws = TRUE tries S3 GeoParquet first; falls back to API on failure.
download_one_taxon_obis <- function(taxon, cfg = config, taxa_cfg = taxa_config,
                                    include_flagged = TRUE, max_age_days = Inf,
                                    force = config$force, use_aws = TRUE, batch_size = 100) {
  tcfg <- taxa_cfg[[taxon]]
  if (is.null(tcfg)) stop(glue("Taxon '{taxon}' not in taxa_config."))
  if (is.null(tcfg$obis_dataset_path)) {
    stop(glue("Taxon '{taxon}' has no obis_dataset_path in taxa_config."))
  }

  crosswalk <- load_crosswalk(taxon, cfg)
  obis_names <- get_obis_names(crosswalk, include_flagged = include_flagged)
  if (length(obis_names) == 0) stop(glue("No species names resolved for '{taxon}'."))

  out_file <- tcfg$obis_dataset_path
  out_dir <- dirname(out_file)
  req_hash <- obis_request_hash(obis_names, cfg)

  summary_row <- function(status, n_records, source = "api") {
    tibble::tibble(
      taxon = taxon, status = status,
      n_species = nrow(crosswalk), n_names = length(obis_names),
      n_records = n_records, output_file = out_file,
      source = source, run_time_utc = as.character(Sys.time())
    )
  }

  if (!isTRUE(force)) {
    manifest <- read_obis_manifest(out_dir)
    reason <- obis_download_stale_reason(manifest, req_hash, out_file, max_age_days)
    if (is.null(reason)) {
      message(glue("[obis] {taxon}: using CACHED download ",
                   "(downloaded {manifest$downloaded_utc} UTC; nothing re-fetched)."))
      return(summary_row("skipped", manifest$n_records, manifest$source %||% "api"))
    }
    message(glue("[obis] {taxon}: downloading ({reason})."))
  }

  wkt <- obis_region_wkt(cfg)
  bbox <- region_bbox(cfg)
  
  # Try AWS S3 GeoParquet first if enabled
  recs <- NULL
  source <- "api"
  if (isTRUE(use_aws)) {
    message(glue("[obis] {taxon}: trying AWS S3 GeoParquet..."))
    recs <- fetch_obis_aws(obis_names, range(cfg$years), bbox, cfg)
    if (!is.null(recs) && nrow(recs) > 0) {
      source <- "aws_geoparquet"
      message(glue("[obis] {taxon}: {nrow(recs)} records from AWS S3 GeoParquet"))
    }
  }
  
  # Fall back to API if S3 failed or returned nothing
  if (is.null(recs) || nrow(recs) == 0) {
    message(glue("[obis] {taxon}: querying API (batch_size = {batch_size})..."))
    source <- "api"
    batches <- split(obis_names, ceiling(seq_along(obis_names) / batch_size))
    recs <- purrr::imap(batches, function(b, i) {
      message(glue("[obis] {taxon}: batch {i}/{length(batches)} ({length(b)} names)"))
      fetch_obis_batch(b, range(cfg$years), wkt)
    }) |>
      dplyr::bind_rows()
  }
  
  recs <- standardize_obis(recs, cfg)

  ensure_dir(out_dir)
  arrow::write_parquet(recs, out_file)
  jsonlite::write_json(
    list(taxon = taxon, request_hash = req_hash, file = basename(out_file),
         n_names = length(obis_names), n_records = nrow(recs), source = source,
         downloaded_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S")),
    obis_manifest_path(out_dir), auto_unbox = TRUE, pretty = TRUE
  )

  summary_row("ok", nrow(recs), source)
}

### Driver over taxa, with per-taxon error capture so one failure doesn't abort.
### Taxa without an `obis_dataset_path` are skipped by default (they are
### GBIF-only); pass them explicitly via `taxa` to get the error instead.
### max_age_days: re-download if the stored file is older than this
### (default: config$obis_max_age_days, else 180; Inf = never).
### use_aws: if TRUE (default), tries AWS S3 GeoParquet first before falling back to API.
run_obis_downloads <- function(taxa = Filter(function(t) !is.null(taxa_cfg[[t]]$obis_dataset_path),
                                             cfg$taxa),
                               cfg = config, taxa_cfg = taxa_config,
                               include_flagged = TRUE,
                               max_age_days = cfg$obis_max_age_days %||% 180,
                               force = config$force, use_aws = TRUE) {
  if (!source_enabled("obis", cfg)) {
    message("[obis] disabled (config$obis_filters$enabled = FALSE); skipping.")
    return(tibble::tibble())
  }
  purrr::map_dfr(taxa, function(taxon) {
    tryCatch(
      download_one_taxon_obis(taxon, cfg, taxa_cfg, include_flagged,
                              max_age_days, force, use_aws),
      error = function(e) tibble::tibble(
        taxon = taxon, status = "error",
        n_species = NA_integer_, n_names = NA_integer_, n_records = NA_integer_,
        output_file = NA_character_, source = NA_character_,
        run_time_utc = as.character(Sys.time()),
        error_message = conditionMessage(e)
      )
    )
  })
}
