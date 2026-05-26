# Extract structured trajectory outputs from a Monocle3 cell_data_set.

ta_extract_trajectory_outputs <- function(
  cds,
  pseudotime,
  pseudotime_bin,
  root_pr_nodes,
  branch_path = NULL,
  output_dir,
  cell_type_col = NULL,
  condition_col = NULL,
  batch_col = NULL,
  donor_col = NULL,
  time_col = NULL
) {
  monocle_dir <- file.path(output_dir, "monocle_outputs")
  graph <- ta_get_pr_graph(cds)
  closest <- ta_get_closest_pr_vertices(cds)
  partitions <- ta_get_partitions(cds)
  clusters <- ta_get_clusters(cds)
  cell_ids <- names(partitions)
  if (is.null(cell_ids) || length(cell_ids) == 0) {
    cell_ids <- colnames(cds)
  }

  degrees <- igraph::degree(graph)
  node_ids <- igraph::V(graph)$name
  node_type <- ifelse(
    node_ids %in% root_pr_nodes,
    "root",
    ifelse(degrees > 2, "branch", ifelse(degrees == 1, "leaf", "internal"))
  )
  nodes_df <- data.frame(
    node_id = node_ids,
    degree = as.integer(degrees),
    node_type = node_type,
    stringsAsFactors = FALSE
  )
  edges_df <- igraph::as_data_frame(graph, what = "edges")
  if (NROW(edges_df) == 0) {
    edges_df <- ta_empty_df(c("source", "target"))
  } else {
    colnames(edges_df)[1:2] <- c("source", "target")
    edges_df <- edges_df[, c("source", "target"), drop = FALSE]
  }

  umap <- SingleCellExperiment::reducedDims(cds)$UMAP
  umap_df <- data.frame(
    cell_id = rownames(umap),
    UMAP_1 = as.numeric(umap[, 1]),
    UMAP_2 = as.numeric(umap[, 2]),
    stringsAsFactors = FALSE
  )

  col_data <- ta_col_data(cds)
  metadata_check <- ta_check_requested_metadata_columns(
    col_data,
    c(cell_type_col, condition_col, batch_col, donor_col, time_col)
  )
  metadata_cols <- metadata_check$present
  branch_path <- branch_path %||% stats::setNames(rep(NA_character_, length(cell_ids)), cell_ids)

  cell_df <- data.frame(
    cell_id = cell_ids,
    pseudotime = as.numeric(pseudotime[cell_ids]),
    pseudotime_bin = as.character(pseudotime_bin[cell_ids]),
    partition = as.character(partitions[cell_ids]),
    cluster = as.character(clusters[cell_ids]),
    closest_pr_node = as.character(closest[cell_ids]),
    branch_path = as.character(branch_path[cell_ids]),
    stringsAsFactors = FALSE
  )
  if (length(metadata_cols) > 0) {
    cell_df <- cbind(cell_df, col_data[cell_ids, metadata_cols, drop = FALSE])
  }

  ta_write_tsv(cell_df, file.path(monocle_dir, "cell_trajectory.tsv"))
  ta_write_tsv(cell_df, file.path(monocle_dir, "cell_pr_node.tsv"))
  ta_write_tsv(umap_df, file.path(monocle_dir, "umap.tsv"))
  ta_write_tsv(nodes_df, file.path(monocle_dir, "principal_graph_nodes.tsv"))
  ta_write_tsv(edges_df, file.path(monocle_dir, "principal_graph_edges.tsv"))

  cells_per_partition <- ta_named_count_list(partitions)
  cells_per_cluster <- ta_named_count_list(clusters)
  cluster_partition_table <- as.data.frame.matrix(table(cluster = clusters, partition = partitions))
  cluster_partition_table <- data.frame(cluster = rownames(cluster_partition_table), cluster_partition_table, row.names = NULL)
  ta_write_tsv(cluster_partition_table, file.path(monocle_dir, "cluster_partition_table.tsv"))

  list(
    cell_trajectory = cell_df,
    umap = umap_df,
    principal_graph_nodes = nodes_df,
    principal_graph_edges = edges_df,
    graph = graph,
    closest_pr_node = closest,
    partitions = partitions,
    clusters = clusters,
    metadata_warnings = metadata_check$warnings,
    missing_metadata_columns = metadata_check$missing,
    trajectory_structure = list(
      num_cells = as.integer(ncol(cds)),
      num_genes = as.integer(nrow(cds)),
      num_clusters = as.integer(length(unique(clusters))),
      num_partitions = as.integer(length(unique(partitions))),
      cells_per_partition = cells_per_partition,
      cells_per_cluster = cells_per_cluster,
      num_branch_points = as.integer(sum(degrees > 2)),
      num_leaf_nodes = as.integer(sum(degrees == 1)),
      branch_point_ids = as.character(node_ids[degrees > 2]),
      leaf_node_ids = as.character(node_ids[degrees == 1])
    )
  )
}
