### Experimental incremental variant of the greedy loop (verification only).
### Sourced after the baseline; redefines greedy_setcover with incremental scoring.
### Keeps `edges` semantics identical to the baseline: batch cells are NOT removed
### from edges (they may be re-chosen), only satisfied species are pruned.

greedy_setcover <- function(year, taxon, cfg = config, params = sampling_params,
                            force = config$force,
                            progress = c("quiet", "terse", "verbose")) {
  progress <- match.arg(progress)
  year <- as.character(year)
  taxa <- resolve_greedy_taxa(taxon, cfg)
  samp <- make_sampling_config(year, cfg = cfg)
  cons_balance_weight <- params$cons_balance_weight
  batch_fraction      <- params$batch_fraction
  deficit_percentile  <- params$deficit_percentile
  output_dir <- build_greedy_output_dir(year, taxon, cfg, taxa)

  updated_ranges <- rbindlist(lapply(taxa, function(t) {
    dt <- as.data.table(read_parquet(samp[[t]]$undersamp_cells_out)); dt[, taxon := t]; dt
  }))
  undersamp_species <- rbindlist(lapply(taxa, function(t) {
    df <- read_csv(samp[[t]]$undersamp_spp_out, show_col_types = FALSE); df$taxon <- t; df
  }), fill = TRUE)
  range_cells_long <- rbindlist(lapply(taxa, function(t) {
    rc <- as.data.table(read_parquet(samp[[t]]$range_cells_path))
    melt(rc, id.vars = "h3_index", variable.name = "species",
         value.name = "in_range")[in_range == 1]
  }))
  CA_hex_grid <- .hex_grid_with_cons(cfg)
  cons_lookup <- unique(as.data.table(sf::st_drop_geometry(CA_hex_grid))[, .(h3_index, cons_area)])

  edges <- as.data.table(updated_ranges)[pres_abs == 1, .(h3_index, species)]
  edges <- merge(edges, cons_lookup, by = "h3_index", all.x = TRUE)
  edges[is.na(cons_area), cons_area := 0L]
  zm_long <- merge(range_cells_long, cons_lookup, by = "h3_index", all.x = TRUE)
  zm_long[is.na(cons_area), cons_area := 0L]
  cons_target <- zm_long[, .(cons_frac_range = sum(cons_area) / .N), by = species]
  spp_dt <- as.data.table(undersamp_species)[, .(species, taxon, cells_needed)]
  spp_dt <- merge(spp_dt, cons_target, by = "species", all.x = TRUE)
  spp_dt[is.na(cons_frac_range), cons_frac_range := 0]

  remaining <- copy(spp_dt); setkey(remaining, species)
  selected <- vector("list", 5000L)
  n_batches <- 0L; n_selected <- 0L
  total_cons_selected <- 0L; total_noncons_selected <- 0L
  setkey(edges, species)
  edges <- edges[species %in% remaining$species]

  ### Per-cell score built once; decremented by removed edges each batch.
  cell_scores <- edges[, .(n_spp = .N, cons_area = cons_area[1]), by = h3_index]

  while (nrow(remaining[cells_needed > 0]) > 0) {
    cell_scores <- cell_scores[n_spp > 0]
    if (nrow(cell_scores) == 0) break

    median_target <- median(remaining[cells_needed > 0, cons_frac_range])
    current_cons_frac <- if ((total_cons_selected + total_noncons_selected) == 0) 0.5 else
      total_cons_selected / (total_cons_selected + total_noncons_selected)
    cons_gap <- median_target - current_cons_frac
    cell_scores[, adj_score := n_spp + cons_balance_weight * n_spp *
                  fifelse(cons_area == 1, cons_gap, -cons_gap)]

    deficit_p <- remaining[cells_needed > 0, quantile(cells_needed, probs = deficit_percentile)]
    batch_size <- max(1L, floor(deficit_p * batch_fraction))
    batch_size <- min(batch_size, nrow(cell_scores))

    setorder(cell_scores, -adj_score)
    batch <- cell_scores[1:batch_size]
    batch_cells <- batch$h3_index

    ### Batch edges for still-short species only (what the baseline counts as hits).
    active_spp <- remaining[cells_needed > 0, species]
    batch_edges <- edges[J(active_spp), nomatch = 0L][h3_index %in% batch_cells]
    spp_hits <- batch_edges[, .(hits = .N), by = species]
    remaining[spp_hits, cells_needed := cells_needed - i.hits, on = "species"]

    batch_cons <- sum(batch$cons_area == 1L)
    total_cons_selected <- total_cons_selected + batch_cons
    total_noncons_selected <- total_noncons_selected + (batch_size - batch_cons)

    n_batches <- n_batches + 1L
    selected[[n_batches]] <- data.table(
      order = seq(n_selected + 1L, n_selected + batch_size),
      h3_index = batch$h3_index, cons_area = batch$cons_area,
      n_spp_covered = batch$n_spp, batch_id = n_batches,
      spp_remaining = nrow(remaining[cells_needed > 0])
    )
    n_selected <- n_selected + batch_size
    .progress_msg(progress, "verbose", glue("    batch {n_batches}: {n_selected} cells selected, ",
                  "{nrow(remaining[cells_needed > 0])} species still short"))

    ### Score update: the baseline recounts edges of active species each batch. Those
    ### counts change only by (a) batch-cell edges of active species that were just
    ### matched, no: the baseline does NOT remove them from edges, so counts only
    ### drop when a species becomes satisfied (its edges leave `edges`).
    satisfied <- remaining[cells_needed <= 0, species]
    if (length(satisfied) > 0) {
      gone <- edges[species %in% satisfied, .(h3_index)]
      drop_counts <- gone[, .(drop = .N), by = h3_index]
      cell_scores[drop_counts, n_spp := n_spp - i.drop, on = "h3_index"]
      edges <- edges[!species %in% satisfied]
    }
  }

  selected_cells <- rbindlist(selected[seq_len(n_batches)])
  species_summary <- merge(spp_dt, remaining[, .(species, final_deficit = cells_needed)],
                           by = "species")
  species_summary[, cells_selected := cells_needed - final_deficit]
  species_summary[, target_met := final_deficit <= 0]
  taxon_summary <- species_summary[, .(n_species = .N, n_at_target = sum(target_met),
    pct_at_target = round(100 * sum(target_met) / .N, 1)), by = taxon]
  list(status = "written", selected_cells = selected_cells, species_summary = species_summary,
       taxon_summary = taxon_summary, output_dir = output_dir)
}
