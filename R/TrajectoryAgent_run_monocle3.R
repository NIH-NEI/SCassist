# Native Monocle3 execution helpers for TrajectoryAgent.

ta_normalize_cds_metadata <- function(cds, input_type) {
  warnings <- list()
  cell_data <- ta_col_data(cds)
  if (is.null(rownames(cell_data)) || !identical(rownames(cell_data), colnames(cds))) {
    if (NROW(cell_data) == ncol(cds)) {
      rownames(cell_data) <- colnames(cds)
      SummarizedExperiment::colData(cds) <- S4Vectors::DataFrame(cell_data)
    } else {
      stop("cell metadata must match colnames(cds).", call. = FALSE)
    }
  }
  row_data <- ta_row_data(cds)
  if (!"gene_short_name" %in% colnames(row_data) ||
    all(is.na(row_data$gene_short_name) | !nzchar(as.character(row_data$gene_short_name)))) {
    row_data$gene_short_name <- rownames(cds)
    rownames(row_data) <- rownames(cds)
    SummarizedExperiment::rowData(cds) <- S4Vectors::DataFrame(row_data)
    warnings <- ta_append_warning(
      warnings,
      "gene_short_name_missing_created_from_gene_ids",
      "gene_short_name was absent or empty; gene IDs were used."
    )
  }
  attr(cds, "trajectory_agent_warnings") <- warnings
  attr(cds, "trajectory_agent_input_type") <- input_type
  cds
}

ta_prepare_cds <- function(
  cds = NULL,
  expression_matrix = NULL,
  cell_metadata = NULL,
  gene_metadata = NULL,
  seurat_object = NULL
) {
  if (!is.null(cds)) {
    if (!inherits(cds, "cell_data_set")) {
      stop("cds must be a Monocle3 cell_data_set.", call. = FALSE)
    }
    return(ta_normalize_cds_metadata(cds, "cds"))
  }

  if (!is.null(seurat_object)) {
    if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
      stop("Seurat input requires SeuratWrappers::as.cell_data_set(). Provide cds or expression_matrix instead.", call. = FALSE)
    }
    return(ta_normalize_cds_metadata(SeuratWrappers::as.cell_data_set(seurat_object), "seurat_object"))
  }

  if (is.null(expression_matrix) || is.null(cell_metadata) || is.null(gene_metadata)) {
    stop("Provide cds, seurat_object, or expression_matrix with cell_metadata and gene_metadata.", call. = FALSE)
  }
  if (is.null(rownames(expression_matrix)) || is.null(colnames(expression_matrix))) {
    stop("expression_matrix must have gene row names and cell column names.", call. = FALSE)
  }
  cell_metadata <- as.data.frame(cell_metadata, stringsAsFactors = FALSE)
  gene_metadata <- as.data.frame(gene_metadata, stringsAsFactors = FALSE)
  if (is.null(rownames(cell_metadata)) || is.null(rownames(gene_metadata))) {
    stop("cell_metadata and gene_metadata must have row names.", call. = FALSE)
  }
  if (!identical(colnames(expression_matrix), rownames(cell_metadata))) {
    stop("Expression matrix columns must match cell_metadata row names.", call. = FALSE)
  }
  if (!identical(rownames(expression_matrix), rownames(gene_metadata))) {
    stop("Expression matrix rows must match gene_metadata row names.", call. = FALSE)
  }
  if (!"gene_short_name" %in% colnames(gene_metadata)) {
    gene_metadata$gene_short_name <- rownames(gene_metadata)
    warnings <- list(ta_warning(
      "gene_short_name_missing_created_from_gene_ids",
      "gene_metadata did not contain gene_short_name; gene IDs were used."
    ))
  } else {
    warnings <- list()
  }
  cds <- monocle3::new_cell_data_set(
    expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )
  attr(cds, "trajectory_agent_warnings") <- warnings
  attr(cds, "trajectory_agent_input_type") <- "expression_matrix"
  cds
}

ta_run_monocle3_workflow <- function(
  cds,
  num_dim = 50,
  alignment_group = NULL,
  residual_model_formula_str = NULL,
  seed = 1234
) {
  set.seed(seed)
  cds <- monocle3::preprocess_cds(cds, num_dim = num_dim)
  warnings <- list()
  align_args <- list(cds = cds)
  col_data <- ta_col_data(cds)
  if (!ta_null_arg(alignment_group) && alignment_group %in% colnames(col_data)) {
    align_args$alignment_group <- alignment_group
  } else if (!ta_null_arg(alignment_group)) {
    warnings <- ta_append_warning(
      warnings,
      paste0("alignment_group_missing:", alignment_group),
      paste0("alignment_group '", alignment_group, "' is absent; align_cds will skip alignment_group."),
      list(column = alignment_group)
    )
  }
  if (!ta_null_arg(residual_model_formula_str)) {
    missing_formula_cols <- ta_formula_missing_columns(residual_model_formula_str, col_data)
    if (length(missing_formula_cols) > 0) {
      warnings <- ta_append_warning(
        warnings,
        paste0("residual_model_formula_missing_columns:", paste(missing_formula_cols, collapse = ",")),
        "residual_model_formula_str references absent metadata column(s); residual model will be skipped.",
        list(columns = missing_formula_cols)
      )
    } else {
      align_args$residual_model_formula_str <- residual_model_formula_str
    }
  }
  if (length(align_args) > 1) {
    cds <- do.call(monocle3::align_cds, align_args)
  }
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
  cds <- monocle3::cluster_cells(cds)
  cds <- monocle3::learn_graph(cds)
  attr(cds, "trajectory_agent_warnings") <- c(attr(cds, "trajectory_agent_warnings") %||% list(), warnings)
  cds
}

ta_order_cells_if_rooted <- function(cds, root_selection_evidence) {
  if (length(root_selection_evidence$root_pr_nodes) == 0) {
    return(list(cds = cds, order_cells_ran = FALSE, warnings = list(
      ta_warning("no_defensible_root_selected", "No defensible root was selected; order_cells was not run.")
    )))
  }
  warnings <- list()
  ordered <- tryCatch({
    monocle3::order_cells(cds, root_pr_nodes = root_selection_evidence$root_pr_nodes)
  }, error = function(e) {
    warnings <<- ta_append_warning(warnings, "order_cells_failed", paste("order_cells failed:", conditionMessage(e)))
    cds
  })
  list(cds = ordered, order_cells_ran = length(warnings) == 0, warnings = warnings)
}
