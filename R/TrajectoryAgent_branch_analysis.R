# Non-interactive branch extraction from the Monocle3 principal graph.

ta_analyze_branches <- function(
  graph,
  closest_pr_node,
  cell_trajectory,
  output_dir = NULL,
  cell_type_col = NULL,
  condition_col = NULL,
  batch_col = NULL,
  donor_col = NULL,
  min_cells_per_branch_path = 20,
  branch_associated_genes = list(),
  branch_associated_modules = list()
) {
  degrees <- igraph::degree(graph)
  node_ids <- igraph::V(graph)$name
  branch_nodes <- node_ids[degrees > 2]
  leaf_nodes <- node_ids[degrees == 1]
  branch_summaries <- list()
  path_rows <- list()
  branch_path_by_cell <- stats::setNames(rep(NA_character_, nrow(cell_trajectory)), cell_trajectory$cell_id)
  flags <- character()

  if (length(branch_nodes) == 0) {
    flags <- c(flags, "no_branch_points")
    branch_paths <- ta_empty_df(c("branch_node_id", "path_id", "cell_count", "finite_pseudotime_cell_count", "reachable_leaf_nodes", "confidence_flags"))
    if (!is.null(output_dir)) {
      ta_write_json(branch_summaries, file.path(output_dir, "monocle_outputs", "branch_summary.json"))
      ta_write_tsv(branch_paths, file.path(output_dir, "monocle_outputs", "branch_paths.tsv"))
    }
    return(list(
      branch_summary = branch_summaries,
      branch_paths = branch_paths,
      branch_path_by_cell = branch_path_by_cell,
      confidence_flags = flags
    ))
  }

  for (branch_node in branch_nodes) {
    neighbor_names <- igraph::V(graph)$name[igraph::neighbors(graph, branch_node)]
    graph_without_branch <- igraph::delete_vertices(graph, branch_node)
    components <- igraph::components(graph_without_branch)$membership
    path_summaries <- list()
    cells_near_branch <- names(closest_pr_node)[closest_pr_node == branch_node]

    for (idx in seq_along(neighbor_names)) {
      neighbor <- neighbor_names[[idx]]
      if (!neighbor %in% names(components)) {
        next
      }
      component_nodes <- names(components)[components == components[[neighbor]]]
      path_id <- paste0(branch_node, "__path_", idx)
      path_cells <- names(closest_pr_node)[closest_pr_node %in% component_nodes]
      branch_path_by_cell[path_cells] <- path_id
      path_df <- cell_trajectory[match(path_cells, cell_trajectory$cell_id), , drop = FALSE]
      finite_count <- sum(is.finite(suppressWarnings(as.numeric(path_df$pseudotime))), na.rm = TRUE)
      reachable_leaf_nodes <- intersect(leaf_nodes, component_nodes)

      dominant_cell_types <- if (!ta_null_arg(cell_type_col) && cell_type_col %in% colnames(path_df)) {
        ta_composition_summary(path_df[[cell_type_col]])
      } else list()
      dominant_conditions <- if (!ta_null_arg(condition_col) && condition_col %in% colnames(path_df)) {
        ta_composition_summary(path_df[[condition_col]])
      } else list()
      dominant_batches <- if (!ta_null_arg(batch_col) && batch_col %in% colnames(path_df)) {
        ta_composition_summary(path_df[[batch_col]])
      } else list()
      dominant_donors <- if (!ta_null_arg(donor_col) && donor_col %in% colnames(path_df)) {
        ta_composition_summary(path_df[[donor_col]])
      } else list()

      terminal_df <- path_df[path_df$closest_pr_node %in% reachable_leaf_nodes, , drop = FALSE]
      finite_pt <- suppressWarnings(as.numeric(path_df$pseudotime))
      if (nrow(terminal_df) == 0 && any(is.finite(finite_pt))) {
        cutoff <- stats::quantile(finite_pt[is.finite(finite_pt)], probs = 0.80, na.rm = TRUE)
        terminal_df <- path_df[is.finite(finite_pt) & finite_pt >= cutoff, , drop = FALSE]
      }
      terminal_annotations <- character()
      if (nrow(terminal_df) > 0 && !ta_null_arg(cell_type_col) && cell_type_col %in% colnames(terminal_df)) {
        terminal_comp <- ta_composition_summary(terminal_df[[cell_type_col]], top_n = 2)
        if (length(terminal_comp) > 0) {
          terminal_annotations <- c(
            terminal_annotations,
            paste0("Terminal cells may be enriched for ", terminal_comp[[1]]$value, " on ", path_id, ".")
          )
        }
      }

      path_flags <- character()
      if (length(path_cells) < min_cells_per_branch_path) {
        path_flags <- c(path_flags, "small_branch_path", "low_cell_count_branch_path")
        flags <- c(flags, "branch_paths_below_min_cell_count")
      }
      path_summaries[[length(path_summaries) + 1]] <- list(
        path_id = path_id,
        cell_count = as.integer(length(path_cells)),
        finite_pseudotime_cell_count = as.integer(finite_count),
        reachable_leaf_nodes = as.character(reachable_leaf_nodes),
        dominant_cell_types = dominant_cell_types,
        dominant_conditions = dominant_conditions,
        terminal_state_annotations = terminal_annotations,
        confidence_flags = unique(path_flags)
      )
      path_rows[[length(path_rows) + 1]] <- data.frame(
        branch_node_id = branch_node,
        path_id = path_id,
        neighbor_node_id = neighbor,
        cell_count = length(path_cells),
        finite_pseudotime_cell_count = finite_count,
        reachable_leaf_nodes = paste(reachable_leaf_nodes, collapse = ","),
        dominant_cell_types = paste(vapply(dominant_cell_types, function(x) x$value, character(1)), collapse = ","),
        dominant_conditions = paste(vapply(dominant_conditions, function(x) x$value, character(1)), collapse = ","),
        dominant_batches = paste(vapply(dominant_batches, function(x) x$value, character(1)), collapse = ","),
        dominant_donors = paste(vapply(dominant_donors, function(x) x$value, character(1)), collapse = ","),
        terminal_state_annotations = paste(terminal_annotations, collapse = "|"),
        confidence_flags = paste(path_flags, collapse = ","),
        stringsAsFactors = FALSE
      )
    }

    supported <- vapply(path_summaries, function(path) path$cell_count >= min_cells_per_branch_path, logical(1))
    adequate_paths <- path_summaries[supported]
    enough_paths <- length(adequate_paths) >= 2
    has_terminal <- any(vapply(path_summaries, function(path) length(path$terminal_state_annotations) > 0, logical(1)))
    dominant_value <- function(path, field) {
      values <- path[[field]]
      if (length(values) == 0) {
        return(NA_character_)
      }
      if (is.list(values) && length(values[[1]]$value) > 0) {
        return(as.character(values[[1]]$value))
      }
      NA_character_
    }
    adequate_cell_types <- stats::na.omit(vapply(adequate_paths, dominant_value, character(1), field = "dominant_cell_types"))
    adequate_conditions <- stats::na.omit(vapply(adequate_paths, dominant_value, character(1), field = "dominant_conditions"))
    has_distinguishable_cell_types <- !ta_null_arg(cell_type_col) && length(unique(adequate_cell_types)) >= 2
    has_distinguishable_conditions <- !ta_null_arg(condition_col) && length(unique(adequate_conditions)) >= 2
    branch_genes <- if (!is.null(branch_associated_genes[[branch_node]])) branch_associated_genes[[branch_node]] else branch_associated_genes
    branch_modules <- if (!is.null(branch_associated_modules[[branch_node]])) branch_associated_modules[[branch_node]] else branch_associated_modules
    has_branch_gene_support <- length(branch_genes) > 0
    has_branch_module_support <- length(branch_modules) > 0
    has_gene_or_module_support <- has_branch_gene_support || has_branch_module_support
    has_any_distinguishability <- has_distinguishable_cell_types || has_distinguishable_conditions || has_gene_or_module_support
    homogeneous_cell_types <- length(adequate_cell_types) >= 2 && length(unique(adequate_cell_types)) == 1
    homogeneous_conditions <- length(adequate_conditions) >= 2 && length(unique(adequate_conditions)) == 1
    branch_flags <- character()
    if (!enough_paths) {
      branch_flags <- c(branch_flags, "no_adequate_branch_paths", "too_few_supported_branch_paths")
    }
    if (!has_terminal) {
      branch_flags <- c(branch_flags, "missing_terminal_annotations", "terminal_states_not_annotated")
      flags <- c(flags, "no_interpretable_terminal_states")
    }
    if (!has_any_distinguishability) {
      branch_flags <- c(branch_flags, "branch_paths_not_distinguishable")
    }
    if (homogeneous_cell_types && (homogeneous_conditions || ta_null_arg(condition_col))) {
      branch_flags <- c(branch_flags, "homogeneous_branch_paths")
    }
    if (!has_gene_or_module_support) {
      branch_flags <- c(branch_flags, "no_branch_associated_genes_or_modules")
    }
    confidence <- "low"
    if (enough_paths && has_any_distinguishability && has_terminal && has_gene_or_module_support) {
      confidence <- "high"
    } else if (enough_paths && has_any_distinguishability && has_terminal) {
      confidence <- "medium"
    }
    flags <- c(flags, branch_flags)
    branch_summaries[[length(branch_summaries) + 1]] <- list(
      branch_node_id = branch_node,
      num_outgoing_paths = as.integer(length(path_summaries)),
      cells_near_branch = as.integer(length(cells_near_branch)),
      path_summaries = path_summaries,
      branch_associated_genes = branch_genes,
      branch_associated_modules = branch_modules,
      confidence = confidence,
      confidence_flags = unique(branch_flags)
    )
  }

  branch_paths <- if (length(path_rows) > 0) do.call(rbind, path_rows) else ta_empty_df(c("branch_node_id", "path_id", "cell_count", "finite_pseudotime_cell_count", "reachable_leaf_nodes", "confidence_flags"))
  if (!is.null(output_dir)) {
    ta_write_json(branch_summaries, file.path(output_dir, "monocle_outputs", "branch_summary.json"))
    ta_write_tsv(branch_paths, file.path(output_dir, "monocle_outputs", "branch_paths.tsv"))
  }
  list(
    branch_summary = branch_summaries,
    branch_paths = branch_paths,
    branch_path_by_cell = branch_path_by_cell,
    confidence_flags = unique(flags)
  )
}
