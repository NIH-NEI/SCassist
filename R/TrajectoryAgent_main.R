#' Run the SCAssist TrajectoryAgent using native R Monocle3
#'
#' @description Runs a non-interactive Monocle3 trajectory workflow, extracts
#' structured evidence, computes QC/meta-metrics and confidence scores, and
#' writes LLM-ready JSON and markdown outputs. Raw Monocle3 tables are not sent
#' to the LLM; they are summarized into an evidence package first.
#'
#' @param cds Optional Monocle3 cell_data_set.
#' @param expression_matrix Optional genes x cells expression matrix.
#' @param cell_metadata Optional cell metadata data.frame with row names matching cells.
#' @param gene_metadata Optional gene metadata data.frame with row names matching genes.
#' @param seurat_object Optional Seurat object; requires SeuratWrappers if used.
#' @param output_dir Output directory.
#' @param num_dim Number of dimensions for preprocess_cds.
#' @param alignment_group Optional metadata column for align_cds.
#' @param residual_model_formula_str Optional residual model formula for align_cds.
#' @param cell_type_col Optional cell type metadata column.
#' @param condition_col Optional condition metadata column.
#' @param batch_col Optional batch metadata column.
#' @param donor_col Optional donor metadata column.
#' @param time_col Optional time metadata column.
#' @param earliest_time_value Optional earliest time value for root selection.
#' @param root_cell_type_col Optional root cell type metadata column.
#' @param root_cell_type_value Optional root cell type value.
#' @param root_condition_col Optional root condition metadata column.
#' @param root_condition_value Optional root condition value.
#' @param root_pr_nodes Optional explicit principal graph root node IDs.
#' @param root_cells Optional explicit root cell IDs.
#' @param q_value_threshold graph_test q-value threshold.
#' @param morans_i_min Minimum Moran's I threshold.
#' @param top_gene_limit Maximum genes used for directionality summaries.
#' @param min_cells_per_partition Minimum cells for rooted partitions.
#' @param min_root_candidate_cells Minimum root candidate cells per partition.
#' @param min_cells_per_branch_path Minimum cells for confident branch paths.
#' @param pseudotime_bins Number of pseudotime bins.
#' @param cores Number of cores for graph_test.
#' @param seed Random seed.
#' @param save_cds Whether to save the final cds RDS.
#'
#' @return A named list containing the structured TrajectoryAgent evidence package.
#'
#' @export
run_trajectory_agent <- function(
  cds = NULL,
  expression_matrix = NULL,
  cell_metadata = NULL,
  gene_metadata = NULL,
  seurat_object = NULL,
  output_dir,
  num_dim = 50,
  alignment_group = NULL,
  residual_model_formula_str = NULL,
  cell_type_col = NULL,
  condition_col = NULL,
  batch_col = NULL,
  donor_col = NULL,
  time_col = NULL,
  earliest_time_value = NULL,
  root_cell_type_col = NULL,
  root_cell_type_value = NULL,
  root_condition_col = NULL,
  root_condition_value = NULL,
  root_pr_nodes = NULL,
  root_cells = NULL,
  q_value_threshold = 0.05,
  morans_i_min = 0,
  top_gene_limit = 100,
  min_cells_per_partition = 30,
  min_root_candidate_cells = 5,
  min_cells_per_branch_path = 20,
  pseudotime_bins = 10,
  cores = 4,
  seed = 1234,
  save_cds = TRUE
) {
  ta_require_packages(c(
    "monocle3",
    "Matrix",
    "methods",
    "igraph",
    "jsonlite",
    "dplyr",
    "tibble",
    "SummarizedExperiment",
    "SingleCellExperiment",
    "S4Vectors"
  ))
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  monocle_dir <- file.path(output_dir, "monocle_outputs")
  dir.create(monocle_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)
  warnings <- list()
  confidence_flags <- character()

  run_config <- list(
    cli_parameters = list(
      num_dim = num_dim,
      alignment_group = alignment_group,
      residual_model_formula_str = residual_model_formula_str,
      cell_type_col = cell_type_col,
      condition_col = condition_col,
      batch_col = batch_col,
      donor_col = donor_col,
      time_col = time_col,
      earliest_time_value = earliest_time_value,
      root_cell_type_col = root_cell_type_col,
      root_cell_type_value = root_cell_type_value,
      root_condition_col = root_condition_col,
      root_condition_value = root_condition_value,
      root_pr_nodes = root_pr_nodes,
      root_cells = root_cells,
      q_value_threshold = q_value_threshold,
      morans_i_min = morans_i_min,
      top_gene_limit = top_gene_limit,
      min_cells_per_partition = min_cells_per_partition,
      min_root_candidate_cells = min_root_candidate_cells,
      min_cells_per_branch_path = min_cells_per_branch_path,
      pseudotime_bins = pseudotime_bins,
      cores = cores,
      seed = seed
    ),
    package_versions = list(
      R = R.version.string,
      monocle3 = ta_package_version("monocle3"),
      Matrix = ta_package_version("Matrix"),
      igraph = ta_package_version("igraph"),
      jsonlite = ta_package_version("jsonlite"),
      dplyr = ta_package_version("dplyr"),
      tibble = ta_package_version("tibble")
    )
  )

  cds <- ta_prepare_cds(cds, expression_matrix, cell_metadata, gene_metadata, seurat_object)
  warnings <- c(warnings, attr(cds, "trajectory_agent_warnings") %||% list())
  run_config$input_type <- attr(cds, "trajectory_agent_input_type") %||% "unknown"
  metadata_check <- ta_check_requested_metadata_columns(
    ta_col_data(cds),
    c(cell_type_col, condition_col, batch_col, donor_col, time_col, alignment_group, root_cell_type_col, root_condition_col)
  )
  warnings <- c(warnings, metadata_check$warnings)
  confidence_flags <- c(confidence_flags, vapply(metadata_check$warnings, function(w) w$flag, character(1)))
  if (!ta_null_arg(root_cell_type_value) && !ta_null_arg(root_cell_type_col) && root_cell_type_col %in% metadata_check$missing) {
    warnings <- ta_append_warning(
      warnings,
      paste0("root_selection_metadata_column_missing:", root_cell_type_col),
      "root_cell_type_value was supplied but root_cell_type_col is absent; this root evidence cannot be used.",
      list(column = root_cell_type_col)
    )
    confidence_flags <- c(confidence_flags, paste0("root_selection_metadata_column_missing:", root_cell_type_col))
  } else if (!ta_null_arg(root_cell_type_value) && ta_null_arg(root_cell_type_col)) {
    warnings <- ta_append_warning(
      warnings,
      "root_selection_metadata_column_missing:root_cell_type_col",
      "root_cell_type_value was supplied without root_cell_type_col; this root evidence cannot be used."
    )
    confidence_flags <- c(confidence_flags, "root_selection_metadata_column_missing:root_cell_type_col")
  }
  if (!ta_null_arg(root_condition_value) && !ta_null_arg(root_condition_col) && root_condition_col %in% metadata_check$missing) {
    warnings <- ta_append_warning(
      warnings,
      paste0("root_selection_metadata_column_missing:", root_condition_col),
      "root_condition_value was supplied but root_condition_col is absent; this root evidence cannot be used.",
      list(column = root_condition_col)
    )
    confidence_flags <- c(confidence_flags, paste0("root_selection_metadata_column_missing:", root_condition_col))
  } else if (!ta_null_arg(root_condition_value) && ta_null_arg(root_condition_col)) {
    warnings <- ta_append_warning(
      warnings,
      "root_selection_metadata_column_missing:root_condition_col",
      "root_condition_value was supplied without root_condition_col; this root evidence cannot be used."
    )
    confidence_flags <- c(confidence_flags, "root_selection_metadata_column_missing:root_condition_col")
  }
  cds <- ta_run_monocle3_workflow(cds, num_dim, alignment_group, residual_model_formula_str, seed)
  workflow_warnings <- attr(cds, "trajectory_agent_warnings") %||% list()
  warnings <- c(warnings, workflow_warnings)
  confidence_flags <- c(confidence_flags, vapply(workflow_warnings, function(w) w$flag, character(1)))

  root_selection <- ta_select_roots(
    cds = cds,
    root_pr_nodes = root_pr_nodes,
    root_cells = root_cells,
    time_col = time_col,
    earliest_time_value = earliest_time_value,
    root_cell_type_col = root_cell_type_col,
    root_cell_type_value = root_cell_type_value,
    root_condition_col = root_condition_col,
    root_condition_value = root_condition_value,
    min_cells_per_partition = min_cells_per_partition,
    min_root_candidate_cells = min_root_candidate_cells
  )
  if (length(root_selection$root_pr_nodes) == 0) {
    confidence_flags <- c(confidence_flags, "no_defensible_root_selected")
  }
  if (length(root_selection$partitions_without_roots) > 0) {
    confidence_flags <- c(confidence_flags, "partitions_without_roots")
  }

  ordered <- ta_order_cells_if_rooted(cds, root_selection)
  cds <- ordered$cds
  warnings <- c(warnings, ordered$warnings)
  run_config$order_cells_ran <- ordered$order_cells_ran
  run_config$root_selection_method <- root_selection$method

  cell_ids <- colnames(cds)
  partitions <- ta_get_partitions(cds)
  pseudotime <- if (ordered$order_cells_ran) {
    pt <- monocle3::pseudotime(cds)
    pt[cell_ids]
  } else {
    stats::setNames(rep(NA_real_, length(cell_ids)), cell_ids)
  }
  names(pseudotime) <- cell_ids
  qc_result <- ta_compute_pseudotime_qc(
    pseudotime,
    partitions[cell_ids],
    pseudotime_bins,
    root_selection$partitions_without_roots
  )
  confidence_flags <- c(confidence_flags, qc_result$confidence_flags)

  col_data <- SummarizedExperiment::colData(cds)
  col_data$pseudotime <- as.numeric(pseudotime[cell_ids])
  col_data$pseudotime_bin <- qc_result$pseudotime_bin[cell_ids]
  col_data$partition <- partitions[cell_ids]
  col_data$cluster <- ta_get_clusters(cds)[cell_ids]
  SummarizedExperiment::colData(cds) <- col_data

  extracted <- ta_extract_trajectory_outputs(
    cds,
    pseudotime,
    qc_result$pseudotime_bin,
    root_selection$root_pr_nodes,
    output_dir = output_dir,
    cell_type_col = cell_type_col,
    condition_col = condition_col,
    batch_col = batch_col,
    donor_col = donor_col,
    time_col = time_col
  )

  branch_result <- ta_analyze_branches(
    graph = extracted$graph,
    closest_pr_node = extracted$closest_pr_node,
    cell_trajectory = extracted$cell_trajectory,
    output_dir = output_dir,
    cell_type_col = cell_type_col,
    condition_col = condition_col,
    batch_col = batch_col,
    donor_col = donor_col,
    min_cells_per_branch_path = min_cells_per_branch_path
  )
  confidence_flags <- c(confidence_flags, branch_result$confidence_flags)

  col_data <- SummarizedExperiment::colData(cds)
  col_data$branch_path <- branch_result$branch_path_by_cell[cell_ids]
  SummarizedExperiment::colData(cds) <- col_data

  extracted <- ta_extract_trajectory_outputs(
    cds,
    pseudotime,
    qc_result$pseudotime_bin,
    root_selection$root_pr_nodes,
    branch_path = branch_result$branch_path_by_cell,
    output_dir = output_dir,
    cell_type_col = cell_type_col,
    condition_col = condition_col,
    batch_col = batch_col,
    donor_col = donor_col,
    time_col = time_col
  )

  dynamic_result <- ta_run_dynamic_genes(
    cds,
    pseudotime,
    qc_result$pseudotime_bin,
    output_dir,
    q_value_threshold,
    morans_i_min,
    top_gene_limit,
    cores
  )
  warnings <- c(warnings, dynamic_result$warnings)
  confidence_flags <- c(confidence_flags, dynamic_result$confidence_flags)

  module_result <- ta_run_gene_modules(
    cds,
    dynamic_result$significant_genes,
    dynamic_result$graph_test,
    dynamic_result$dynamic_gene_trends,
    extracted$cell_trajectory,
    output_dir,
    cell_type_col,
    condition_col
  )
  warnings <- c(warnings, module_result$warnings)
  confidence_flags <- c(confidence_flags, module_result$confidence_flags)

  metadata_result <- ta_compute_metadata_associations(
    extracted$cell_trajectory,
    cell_type_col,
    condition_col,
    batch_col,
    donor_col,
    time_col,
    earliest_time_value
  )
  confidence_flags <- c(confidence_flags, metadata_result$confidence_flags)

  confidence <- ta_score_trajectory_confidence(
    root_selection,
    qc_result$pseudotime_qc,
    dynamic_result$dynamic_gene_summary,
    module_result$gene_module_summary,
    metadata_result$metadata_associations,
    branch_result$branch_summary,
    confidence_flags
  )

  evidence <- ta_create_evidence_package(
    run_config = run_config,
    trajectory_structure = extracted$trajectory_structure,
    root_selection_evidence = root_selection,
    pseudotime_qc = qc_result$pseudotime_qc,
    branch_summary = branch_result$branch_summary,
    dynamic_gene_summary = dynamic_result$dynamic_gene_summary,
    gene_module_summary = module_result$gene_module_summary,
    metadata_associations = metadata_result$metadata_associations,
    confidence_flags = confidence$flags,
    confidence_score = confidence$score,
    confidence_label = confidence$label
  )
  evidence$llm_interpretation_input <- ta_deterministic_interpretation(evidence)

  monocle_run_summary <- list(
    trajectory_structure = evidence$trajectory_structure,
    root_selection_evidence = evidence$root_selection_evidence,
    pseudotime_qc = evidence$pseudotime_qc,
    dynamic_gene_summary = evidence$dynamic_gene_summary,
    gene_module_summary = evidence$gene_module_summary,
    confidence_flags = evidence$confidence_flags,
    confidence_score = evidence$confidence_score,
    confidence_label = evidence$confidence_label
  )
  ta_write_json(run_config, file.path(monocle_dir, "run_config.json"))
  ta_write_json(monocle_run_summary, file.path(monocle_dir, "monocle_run_summary.json"))
  ta_write_json(warnings, file.path(monocle_dir, "warnings.json"))
  ta_write_trajectory_reports(evidence, output_dir)
  if (isTRUE(save_cds)) {
    saveRDS(cds, file.path(monocle_dir, "trajectory_agent_cds.rds"))
  }
  evidence
}
