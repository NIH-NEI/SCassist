# Metadata by pseudotime summaries and simple confounding diagnostics.

ta_composition_by_bin <- function(cell_trajectory, bin_col, value_col) {
  if (ta_null_arg(value_col) || !value_col %in% colnames(cell_trajectory) || !bin_col %in% colnames(cell_trajectory)) {
    return(list())
  }
  bins <- as.character(cell_trajectory[[bin_col]])
  values <- as.character(cell_trajectory[[value_col]])
  out <- list()
  for (bin in sort(unique(bins[!is.na(bins) & nzchar(bins)]))) {
    out[[bin]] <- ta_composition_summary(values[bins == bin])
  }
  out
}

ta_compute_metadata_associations <- function(
  cell_trajectory,
  cell_type_col = NULL,
  condition_col = NULL,
  batch_col = NULL,
  donor_col = NULL,
  time_col = NULL,
  earliest_time_value = NULL
) {
  flags <- character()
  confounding <- character()
  associations <- list(
    cell_type_by_pseudotime = list(by_bin = ta_composition_by_bin(cell_trajectory, "pseudotime_bin", cell_type_col)),
    condition_by_pseudotime = list(by_bin = ta_composition_by_bin(cell_trajectory, "pseudotime_bin", condition_col)),
    batch_by_pseudotime = list(by_bin = ta_composition_by_bin(cell_trajectory, "pseudotime_bin", batch_col)),
    donor_by_pseudotime = list(by_bin = ta_composition_by_bin(cell_trajectory, "pseudotime_bin", donor_col)),
    time_by_pseudotime = list(by_bin = ta_composition_by_bin(cell_trajectory, "pseudotime_bin", time_col)),
    possible_batch_or_donor_confounding = character()
  )

  pseudotime <- suppressWarnings(as.numeric(cell_trajectory$pseudotime))
  finite <- is.finite(pseudotime)
  for (entry in list(list(column = batch_col, label = "batch"), list(column = donor_col, label = "donor"))) {
    column <- entry$column
    label <- entry$label
    if (!ta_null_arg(column) && column %in% colnames(cell_trajectory)) {
      comp <- ta_composition_by_bin(cell_trajectory, "pseudotime_bin", column)
      dominant <- 0
      for (bin in names(comp)) {
        if (length(comp[[bin]]) > 0) {
          dominant <- max(dominant, comp[[bin]][[1]]$fraction)
        }
      }
      if (dominant >= 0.85) {
        flag <- paste0(label, "_strongly_associated_with_pseudotime")
        flags <- c(flags, flag)
        confounding <- c(confounding, flag)
      }
      if (sum(finite) >= 10 && length(unique(cell_trajectory[[column]][finite])) > 1) {
        kw <- tryCatch(stats::kruskal.test(pseudotime[finite] ~ as.factor(cell_trajectory[[column]][finite])), error = function(e) NULL)
        target <- if (label == "batch") "batch_by_pseudotime" else "donor_by_pseudotime"
        associations[[target]]$kruskal_p_value <- if (is.null(kw)) NA_real_ else kw$p.value
        if (!is.null(kw) && is.finite(kw$p.value) && kw$p.value < 0.001) {
          flag <- paste0(label, "_pseudotime_kruskal_p_lt_0.001")
          flags <- c(flags, flag)
          confounding <- c(confounding, flag)
        }
      }
    }
  }

  if (!ta_null_arg(time_col) && time_col %in% colnames(cell_trajectory) && sum(finite) >= 5) {
    values <- cell_trajectory[[time_col]][finite]
    numeric_time <- suppressWarnings(as.numeric(values))
    if (sum(!is.na(numeric_time)) < max(3, length(values) / 2)) {
      ordered <- sort(unique(as.character(values)))
      if (!ta_null_arg(earliest_time_value) && as.character(earliest_time_value) %in% ordered) {
        ordered <- c(as.character(earliest_time_value), setdiff(ordered, as.character(earliest_time_value)))
      }
      numeric_time <- match(as.character(values), ordered)
    }
    if (length(unique(numeric_time[!is.na(numeric_time)])) > 1) {
      rho <- suppressWarnings(stats::cor(numeric_time, pseudotime[finite], method = "spearman", use = "complete.obs"))
      associations$time_by_pseudotime$spearman_rho <- as.numeric(rho)
      if (is.na(rho) || rho < 0.25) {
        flags <- c(flags, "pseudotime_direction_conflicts_with_time_metadata")
      }
    }
  }

  associations$possible_batch_or_donor_confounding <- unique(confounding)
  list(metadata_associations = associations, confidence_flags = unique(flags))
}
