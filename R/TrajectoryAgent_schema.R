# Schema construction and validation helpers for TrajectoryAgent evidence.

ta_create_evidence_package <- function(
  run_config,
  trajectory_structure,
  root_selection_evidence,
  pseudotime_qc,
  branch_summary,
  dynamic_gene_summary,
  gene_module_summary,
  metadata_associations,
  confidence_flags,
  confidence_score,
  confidence_label,
  llm_interpretation_input = list()
) {
  list(
    agent_name = "TrajectoryAgent",
    backend_tool = "Monocle3",
    run_config = run_config,
    trajectory_structure = trajectory_structure,
    root_selection_evidence = root_selection_evidence,
    pseudotime_qc = pseudotime_qc,
    branch_summary = branch_summary,
    dynamic_gene_summary = dynamic_gene_summary,
    gene_module_summary = gene_module_summary,
    metadata_associations = metadata_associations,
    confidence_flags = unique(as.character(confidence_flags)),
    confidence_score = as.numeric(confidence_score),
    confidence_label = as.character(confidence_label),
    llm_interpretation_input = llm_interpretation_input
  )
}

ta_validate_evidence_schema <- function(evidence) {
  required <- c(
    "agent_name",
    "backend_tool",
    "run_config",
    "trajectory_structure",
    "root_selection_evidence",
    "pseudotime_qc",
    "branch_summary",
    "dynamic_gene_summary",
    "gene_module_summary",
    "metadata_associations",
    "confidence_flags",
    "confidence_score",
    "confidence_label",
    "llm_interpretation_input"
  )
  missing <- setdiff(required, names(evidence))
  if (length(missing) > 0) {
    stop("Evidence package is missing fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!identical(evidence$agent_name, "TrajectoryAgent")) {
    stop("agent_name must be 'TrajectoryAgent'", call. = FALSE)
  }
  if (!identical(evidence$backend_tool, "Monocle3")) {
    stop("backend_tool must be 'Monocle3'", call. = FALSE)
  }
  invisible(TRUE)
}

ta_empty_dynamic_gene_summary <- function(q_value_threshold = 0.05, morans_i_min = 0) {
  list(
    num_tested_genes = 0L,
    num_significant_genes = 0L,
    q_value_threshold = q_value_threshold,
    morans_i_min = morans_i_min,
    top_genes_by_morans_I = list(),
    top_genes_by_q_value = list(),
    increasing_genes = list(),
    decreasing_genes = list(),
    transient_genes = list()
  )
}

ta_empty_gene_module_summary <- function() {
  list(
    num_modules = 0L,
    modules = list(),
    module_activity_by_cell_type = list(),
    module_activity_by_condition = list(),
    module_activity_by_partition = list(),
    module_activity_by_pseudotime_bin = list(),
    module_activity_by_branch = list()
  )
}

ta_empty_metadata_associations <- function() {
  list(
    cell_type_by_pseudotime = list(),
    condition_by_pseudotime = list(),
    batch_by_pseudotime = list(),
    donor_by_pseudotime = list(),
    time_by_pseudotime = list(),
    possible_batch_or_donor_confounding = character()
  )
}
