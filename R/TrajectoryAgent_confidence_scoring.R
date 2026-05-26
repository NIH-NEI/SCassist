# Deterministic confidence scoring for TrajectoryAgent evidence.

ta_score_trajectory_confidence <- function(
  root_selection_evidence,
  pseudotime_qc,
  dynamic_gene_summary,
  gene_module_summary,
  metadata_associations,
  branch_summary,
  confidence_flags = character()
) {
  score <- 1.0
  flags <- unique(as.character(confidence_flags))

  if (length(root_selection_evidence$root_pr_nodes) == 0) {
    score <- score - 0.30
    flags <- c(flags, "no_defensible_root_selected")
  }
  frac_inf <- pseudotime_qc$fraction_infinite_pseudotime %||% 0
  if (frac_inf > 0.20) {
    score <- score - 0.20
    flags <- c(flags, "high_infinite_pseudotime_fraction")
  } else if (frac_inf > 0.05) {
    score <- score - 0.10
    flags <- c(flags, "moderate_infinite_pseudotime_fraction")
  }
  if (length(root_selection_evidence$partitions_without_roots) > 0) {
    score <- score - 0.15
    flags <- c(flags, "partitions_without_roots")
  }
  if (length(metadata_associations$possible_batch_or_donor_confounding) > 0) {
    score <- score - 0.20
    flags <- c(flags, "batch_or_donor_strongly_associated_with_pseudotime")
  }
  sig <- dynamic_gene_summary$num_significant_genes %||% 0L
  if (sig == 0) {
    score <- score - 0.15
    flags <- c(flags, "no_trajectory_variable_genes")
  } else if (sig < 25) {
    score <- score - 0.10
    flags <- c(flags, "few_trajectory_variable_genes")
  }
  has_terminal <- any(vapply(branch_summary, function(branch) {
    paths <- branch$path_summaries %||% list()
    any(vapply(paths, function(path) length(path$terminal_state_annotations %||% character()) > 0, logical(1)))
  }, logical(1)))
  if (!has_terminal) {
    score <- score - 0.10
    flags <- c(flags, "no_interpretable_terminal_states")
  }
  small_paths <- any(vapply(branch_summary, function(branch) {
    paths <- branch$path_summaries %||% list()
    any(vapply(paths, function(path) {
      path_flags <- path$confidence_flags %||% character()
      "low_cell_count_branch_path" %in% path_flags || "small_branch_path" %in% path_flags
    }, logical(1)))
  }, logical(1)))
  if (small_paths) {
    score <- score - 0.10
    flags <- c(flags, "branch_paths_below_min_cell_count")
  }
  if ("pseudotime_direction_conflicts_with_time_metadata" %in% flags) {
    score <- score - 0.10
  }
  if ((gene_module_summary$num_modules %||% 0L) == 0) {
    score <- score - 0.05
    flags <- c(flags, "gene_modules_unavailable")
  }

  score <- max(0, min(1, score))
  label <- if (score >= 0.75) {
    "high"
  } else if (score >= 0.45) {
    "medium"
  } else {
    "low"
  }
  list(score = as.numeric(round(score, 3)), label = label, flags = unique(flags))
}
