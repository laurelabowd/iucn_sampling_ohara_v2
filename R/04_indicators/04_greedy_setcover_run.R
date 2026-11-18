### =============================================================================
### R/04_indicators/04_greedy_setcover_run.R  —  driver
### =============================================================================
### Runs the greedy set-cover per taxon (and optionally "combined") for a year.
###
### `force = TRUE` recomputes; `progress` = "terse" (default) | "verbose" | "quiet".
###
### Requires config + sampling_config + R/00_utils.R.

source(here::here("R", "04_indicators", "04a_greedy_setcover_fn.R"))

run_greedy_setcover <- function(taxa = config$taxa, year = config$focal_year,
                                cfg = config, params = sampling_params,
                                combined = FALSE, force = config$force,
                                progress = c("terse", "verbose", "quiet")) {
  progress <- match.arg(progress)
  targets <- if (isTRUE(combined)) "combined" else taxa
  res <- lapply(targets, function(t) greedy_setcover(year, t, cfg, params, force, progress))
  invisible(res)
}
