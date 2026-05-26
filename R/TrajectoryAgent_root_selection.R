# Partition-aware root selection for native Monocle3 TrajectoryAgent runs.

ta_get_closest_pr_vertices <- function(cds) {
  if (is.list(cds) && !is.null(cds$closest_pr_vertex)) {
    closest <- as.character(cds$closest_pr_vertex)
    return(closest)
  }
  closest <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  if (is.null(closest)) {
    return(stats::setNames(character(), character()))
  }
  closest_vec <- as.character(closest[, 1])
  names(closest_vec) <- rownames(closest)
  closest_vec
}

ta_get_partitions <- function(cds) {
  if (is.list(cds) && !is.null(cds$partitions)) {
    return(as.character(cds$partitions))
  }
  parts <- as.character(monocle3::partitions(cds, reduction_method = "UMAP"))
  names(parts) <- names(monocle3::partitions(cds, reduction_method = "UMAP"))
  parts
}

ta_get_clusters <- function(cds) {
  if (is.list(cds) && !is.null(cds$clusters)) {
    return(as.character(cds$clusters))
  }
  clusters <- as.character(monocle3::clusters(cds, reduction_method = "UMAP"))
  names(clusters) <- names(monocle3::clusters(cds, reduction_method = "UMAP"))
  clusters
}

ta_get_pr_graph <- function(cds) {
  if (is.list(cds) && !is.null(cds$graph)) {
    return(cds$graph)
  }
  monocle3::principal_graph(cds)[["UMAP"]]
}

ta_get_partition_for_pr_node <- function(cds, pr_node) {
  closest <- ta_get_closest_pr_vertices(cds)
  parts <- ta_get_partitions(cds)
  cells <- intersect(names(closest)[closest == pr_node], names(parts))
  if (length(cells) == 0) {
    return(NA_character_)
  }
  names(sort(table(parts[cells]), decreasing = TRUE))[1]
}

ta_select_roots_from_candidate_cells <- function(
  cds,
  candidate_cell_ids,
  min_cells_per_partition,
  min_root_candidate_cells
) {
  closest <- ta_get_closest_pr_vertices(cds)
  parts <- ta_get_partitions(cds)
  candidate_cell_ids <- intersect(as.character(candidate_cell_ids), names(closest))
  selected_roots <- character()
  root_partitions <- character()
  partitions_without_roots <- character()

  for (partition_id in sort(unique(parts))) {
    partition_cells <- names(parts)[parts == partition_id]
    if (length(partition_cells) < min_cells_per_partition) {
      partitions_without_roots <- c(partitions_without_roots, partition_id)
      next
    }
    candidates <- intersect(candidate_cell_ids, partition_cells)
    if (length(candidates) < min_root_candidate_cells) {
      partitions_without_roots <- c(partitions_without_roots, partition_id)
      next
    }
    node_counts <- sort(table(closest[candidates]), decreasing = TRUE)
    selected_roots <- c(selected_roots, names(node_counts)[1])
    root_partitions <- c(root_partitions, partition_id)
  }

  list(
    root_pr_nodes = unique(selected_roots),
    root_partitions = unique(root_partitions),
    partitions_without_roots = unique(partitions_without_roots),
    candidate_root_cell_count = length(candidate_cell_ids)
  )
}

ta_summarize_root_selection <- function(
  cds,
  selected_root_pr_nodes,
  method,
  metadata_field_used = NULL,
  metadata_value_used = NULL,
  candidate_root_cell_count = 0L,
  min_cells_per_partition = 30
) {
  parts <- ta_get_partitions(cds)
  partition_sizes <- table(parts)
  eligible_partitions <- names(partition_sizes)[partition_sizes >= min_cells_per_partition]
  root_partitions <- unique(stats::na.omit(vapply(
    selected_root_pr_nodes,
    function(node) ta_get_partition_for_pr_node(cds, node),
    character(1)
  )))
  partitions_without_roots <- setdiff(eligible_partitions, root_partitions)
  confidence <- "low"
  notes <- character()
  if (length(selected_root_pr_nodes) > 0 && length(partitions_without_roots) == 0) {
    confidence <- if (method %in% c("explicit_root_pr_nodes", "explicit_root_cells")) "high" else "medium"
  } else if (length(selected_root_pr_nodes) > 0) {
    confidence <- "medium"
    notes <- c(notes, "Some partitions lack selected roots and may receive infinite pseudotime.")
  } else {
    notes <- c(notes, "No defensible root was selected; order_cells was not run.")
  }
  list(
    method = method,
    metadata_field_used = metadata_field_used,
    metadata_value_used = metadata_value_used,
    root_pr_nodes = as.character(selected_root_pr_nodes),
    root_partitions = as.character(root_partitions),
    partitions_without_roots = as.character(partitions_without_roots),
    candidate_root_cell_count = as.integer(candidate_root_cell_count),
    confidence = confidence,
    notes = notes
  )
}

ta_select_roots <- function(
  cds,
  root_pr_nodes = NULL,
  root_cells = NULL,
  time_col = NULL,
  earliest_time_value = NULL,
  root_cell_type_col = NULL,
  root_cell_type_value = NULL,
  root_condition_col = NULL,
  root_condition_value = NULL,
  min_cells_per_partition = 30,
  min_root_candidate_cells = 5
) {
  graph <- ta_get_pr_graph(cds)
  graph_nodes <- igraph::V(graph)$name
  col_data <- ta_col_data(cds)
  selected <- character()
  method <- "none"
  field <- NULL
  value <- NULL
  candidate_count <- 0L

  if (!ta_null_arg(root_pr_nodes)) {
    selected <- intersect(as.character(root_pr_nodes), graph_nodes)
    method <- "explicit_root_pr_nodes"
    candidate_count <- length(selected)
  } else if (!ta_null_arg(root_cells)) {
    result <- ta_select_roots_from_candidate_cells(
      cds,
      root_cells,
      min_cells_per_partition,
      min_root_candidate_cells
    )
    selected <- result$root_pr_nodes
    method <- "explicit_root_cells"
    candidate_count <- result$candidate_root_cell_count
  } else if (!ta_null_arg(time_col) && !ta_null_arg(earliest_time_value) && time_col %in% colnames(col_data)) {
    candidates <- rownames(col_data)[as.character(col_data[[time_col]]) == as.character(earliest_time_value)]
    result <- ta_select_roots_from_candidate_cells(cds, candidates, min_cells_per_partition, min_root_candidate_cells)
    selected <- result$root_pr_nodes
    method <- "earliest_time_point"
    field <- time_col
    value <- earliest_time_value
    candidate_count <- result$candidate_root_cell_count
  } else if (!ta_null_arg(root_cell_type_col) && !ta_null_arg(root_cell_type_value) && root_cell_type_col %in% colnames(col_data)) {
    candidates <- rownames(col_data)[as.character(col_data[[root_cell_type_col]]) == as.character(root_cell_type_value)]
    result <- ta_select_roots_from_candidate_cells(cds, candidates, min_cells_per_partition, min_root_candidate_cells)
    selected <- result$root_pr_nodes
    method <- "starting_cell_type"
    field <- root_cell_type_col
    value <- root_cell_type_value
    candidate_count <- result$candidate_root_cell_count
  } else if (!ta_null_arg(root_condition_col) && !ta_null_arg(root_condition_value) && root_condition_col %in% colnames(col_data)) {
    candidates <- rownames(col_data)[as.character(col_data[[root_condition_col]]) == as.character(root_condition_value)]
    result <- ta_select_roots_from_candidate_cells(cds, candidates, min_cells_per_partition, min_root_candidate_cells)
    selected <- result$root_pr_nodes
    method <- "starting_condition"
    field <- root_condition_col
    value <- root_condition_value
    candidate_count <- result$candidate_root_cell_count
  }

  ta_summarize_root_selection(
    cds = cds,
    selected_root_pr_nodes = unique(selected),
    method = method,
    metadata_field_used = field,
    metadata_value_used = value,
    candidate_root_cell_count = candidate_count,
    min_cells_per_partition = min_cells_per_partition
  )
}
