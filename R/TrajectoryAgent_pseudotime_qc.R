# Pseudotime quality control summaries for TrajectoryAgent.

ta_compute_pseudotime_qc <- function(
  pseudotime,
  partitions,
  pseudotime_bins = 10,
  partitions_without_roots = character()
) {
  pseudotime <- as.numeric(pseudotime)
  finite_mask <- is.finite(pseudotime)
  infinite_mask <- is.infinite(pseudotime)
  bins <- ta_make_pseudotime_bins(pseudotime, pseudotime_bins)
  finite_count <- sum(finite_mask, na.rm = TRUE)
  infinite_count <- sum(infinite_mask, na.rm = TRUE)
  total_cells <- length(pseudotime)
  fraction_infinite <- if (total_cells > 0) infinite_count / total_cells else 0
  pseudotime_range <- if (finite_count > 0) as.numeric(range(pseudotime[finite_mask], na.rm = TRUE)) else NULL
  cells_per_bin <- ta_named_count_list(bins[!is.na(bins)])
  partitions_with_inf <- character()
  if (length(partitions) == length(pseudotime) && infinite_count > 0) {
    partitions_with_inf <- sort(unique(as.character(partitions[infinite_mask])))
  }

  flags <- character()
  if (finite_count == 0) {
    flags <- c(flags, "no_pseudotime_available")
  }
  if (fraction_infinite > 0.20) {
    flags <- c(flags, "high_infinite_pseudotime_fraction")
  } else if (fraction_infinite > 0.05) {
    flags <- c(flags, "moderate_infinite_pseudotime_fraction")
  }

  list(
    pseudotime_qc = list(
      finite_pseudotime_cells = as.integer(finite_count),
      infinite_pseudotime_cells = as.integer(infinite_count),
      fraction_infinite_pseudotime = as.numeric(fraction_infinite),
      pseudotime_range = pseudotime_range,
      cells_per_pseudotime_bin = cells_per_bin,
      partitions_with_infinite_pseudotime = partitions_with_inf,
      partitions_without_roots = as.character(partitions_without_roots)
    ),
    pseudotime_bin = bins,
    confidence_flags = unique(flags)
  )
}
