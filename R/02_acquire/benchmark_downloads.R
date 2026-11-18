### =============================================================================
### R/02_acquire/benchmark_downloads.R  —  compare old vs. new download speeds
### =============================================================================
###
### Benchmarks GBIF batch-submit and OBIS batch-size changes on a small sample.
### Requires config, taxa_config, and both download scripts sourced.
###
### Usage:
###   source(here("R/02_acquire/benchmark_downloads.R"))
###   benchmark_downloads(sample_taxa = "mammal", obis_batch_sizes = c(25, 100))

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(glue)
  library(tibble)
})

source(here("config", "config.R"))
source(here("config", "taxa_config.R"))
source(here("R", "00_utils.R"))
source(here("R", "02_acquire", "01_download_gbif.R"))
source(here("R", "02_acquire", "02_download_obis.R"))

### Benchmark GBIF: current approach
benchmark_gbif <- function(taxa = "mammal", force = FALSE) {
  cat("\n=== GBIF Benchmark (current approach) ===\n")
  
  cat("Timing current sequential approach...\n")
  t_current <- system.time({
    result_current <- run_gbif_downloads(taxa = taxa, force = force, wait = TRUE)
  })
  
  cat("Result:\n")
  print(result_current)
  cat(glue("\nCurrent approach took {round(t_current[3], 1)}s elapsed\n"))
  
  # Return timing
  tibble(
    approach = "current (sequential wait)",
    taxa = paste(taxa, collapse = ", "),
    elapsed_sec = round(t_current[3], 1),
    user_sec = round(t_current[1], 1),
    system_sec = round(t_current[2], 1)
  )
}

### Benchmark OBIS: batch size 100 (new) vs. 25 (simulated old)
benchmark_obis <- function(taxa = NULL, obis_batch_sizes = c(25, 100),
                          force = FALSE, cfg = config, taxa_cfg = taxa_config) {
  if (is.null(taxa)) {
    taxa <- Filter(function(t) !is.null(taxa_cfg[[t]]$obis_dataset_path), cfg$taxa)
  }
  
  cat("\n=== OBIS Benchmark ===\n")
  
  results <- purrr::map_dfr(obis_batch_sizes, function(bs) {
    cat(glue("\nTiming with batch_size = {bs}...\n"))
    t <- system.time({
      result <- purrr::map_dfr(taxa, function(taxon) {
        tryCatch(
          download_one_taxon_obis(taxon, cfg, taxa_cfg,
                                  include_flagged = TRUE, max_age_days = Inf,
                                  force = force, batch_size = bs),
          error = function(e) tibble(
            taxon = taxon, status = "error", n_species = NA_integer_,
            n_names = NA_integer_, n_records = NA_integer_,
            output_file = NA_character_, run_time_utc = as.character(Sys.time()),
            error_message = conditionMessage(e)
          )
        )
      })
    })
    
    cat("Result:\n")
    print(result)
    
    tibble(
      batch_size = bs,
      taxa = paste(taxa, collapse = ", "),
      elapsed_sec = round(t[3], 1),
      user_sec = round(t[1], 1),
      system_sec = round(t[2], 1),
      approx_speedup = paste0("1.0x")  # Will fill in after comparison
    )
  })
  
  # Compute speedup
  if (nrow(results) > 1) {
    baseline <- results$elapsed_sec[1]
    results <- mutate(results,
      approx_speedup = glue("{round(baseline / elapsed_sec, 2)}x")
    )
  }
  
  results
}

### Combined benchmark
benchmark_downloads <- function(gbif_taxa = "mammal", obis_taxa = NULL,
                               obis_batch_sizes = c(25, 100), force = FALSE) {
  cat("\n################################################################################\n")
  cat("# Download Benchmark: GBIF Batch-Submit & OBIS Batch Size\n")
  cat("################################################################################\n")
  
  gbif_res <- if (!is.null(gbif_taxa)) benchmark_gbif(gbif_taxa, force) else tibble()
  obis_res <- if (!is.null(obis_taxa)) benchmark_obis(obis_taxa, obis_batch_sizes, force) else tibble()
  
  cat("\n################################################################################\n")
  cat("# Summary\n")
  cat("################################################################################\n")
  
  if (nrow(gbif_res) > 0) {
    cat("\nGBIF Results:\n")
    print(gbif_res)
  }
  
  if (nrow(obis_res) > 0) {
    cat("\nOBIS Results:\n")
    print(obis_res)
  }
  
  invisible(list(gbif = gbif_res, obis = obis_res))
}
