# Gene module detection and module activity summaries.

ta_activity_to_summary <- function(df) {
  if (is.null(df) || NROW(df) == 0) {
    return(list())
  }
  list(records = ta_records_from_df(df, 200))
}

ta_run_gene_modules <- function(
  cds,
  significant_genes,
  graph_test_df,
  dynamic_gene_trends_df,
  cell_trajectory,
  output_dir,
  cell_type_col = NULL,
  condition_col = NULL
) {
  warnings <- list()
  flags <- character()
  monocle_dir <- file.path(output_dir, "monocle_outputs")
  module_cols <- c("gene_id", "gene_short_name", "module", "q_value", "morans_I", "trend_category")
  activity_files <- c(
    cell_type = "module_activity_by_cell_type.tsv",
    condition = "module_activity_by_condition.tsv",
    partition = "module_activity_by_partition.tsv",
    pseudotime_bin = "module_activity_by_pseudotime_bin.tsv",
    branch = "module_activity_by_branch_path.tsv"
  )

  if (is.null(significant_genes) || NROW(significant_genes) < 10) {
    flags <- c(flags, "too_few_genes_for_gene_modules", "no_gene_modules")
    ta_write_tsv(ta_empty_df(module_cols), file.path(monocle_dir, "gene_modules.tsv"))
    for (file in activity_files) ta_write_tsv(ta_empty_df(c("module")), file.path(monocle_dir, file))
    return(list(
      gene_modules = ta_empty_df(module_cols),
      gene_module_summary = ta_empty_gene_module_summary(),
      warnings = warnings,
      confidence_flags = flags
    ))
  }

  module_df <- tryCatch({
    monocle3::find_gene_modules(cds[significant_genes$gene_id, ], resolution = c(10^seq(-6, -1)))
  }, error = function(e) {
    warnings <<- ta_append_warning(warnings, "find_gene_modules_failed", paste("find_gene_modules failed:", conditionMessage(e)))
    NULL
  })
  if (is.null(module_df) || NROW(module_df) == 0) {
    flags <- c(flags, "no_gene_modules")
    ta_write_tsv(ta_empty_df(module_cols), file.path(monocle_dir, "gene_modules.tsv"))
    for (file in activity_files) ta_write_tsv(ta_empty_df(c("module")), file.path(monocle_dir, file))
    return(list(
      gene_modules = ta_empty_df(module_cols),
      gene_module_summary = ta_empty_gene_module_summary(),
      warnings = warnings,
      confidence_flags = flags
    ))
  }

  module_df <- as.data.frame(module_df, stringsAsFactors = FALSE)
  if (!"gene_id" %in% colnames(module_df)) {
    module_df$gene_id <- if ("id" %in% colnames(module_df)) module_df$id else rownames(module_df)
  }
  row_data <- ta_row_data(cds)
  row_data$gene_id <- rownames(row_data)
  module_df <- ta_normalize_gene_short_name(module_df, row_data, gene_id_col = "gene_id")
  module_df <- merge(module_df, graph_test_df[, intersect(c("gene_id", "q_value", "morans_I"), colnames(graph_test_df)), drop = FALSE], by = "gene_id", all.x = TRUE)
  if (!is.null(dynamic_gene_trends_df) && NROW(dynamic_gene_trends_df) > 0) {
    module_df <- merge(module_df, dynamic_gene_trends_df[, intersect(c("gene_id", "trend_category"), colnames(dynamic_gene_trends_df)), drop = FALSE], by = "gene_id", all.x = TRUE)
  } else {
    module_df$trend_category <- NA_character_
  }
  ta_write_tsv(module_df[, intersect(module_cols, colnames(module_df)), drop = FALSE], file.path(monocle_dir, "gene_modules.tsv"))

  activity <- list()
  write_activity <- function(column, key, file) {
    if (ta_null_arg(column) || !column %in% colnames(cell_trajectory)) {
      ta_write_tsv(ta_empty_df(c("module")), file.path(monocle_dir, file))
      return(list())
    }
    cell_group_df <- data.frame(
      cell = cell_trajectory$cell_id,
      cell_group = as.character(cell_trajectory[[column]]),
      stringsAsFactors = FALSE
    )
    act <- tryCatch(monocle3::aggregate_gene_expression(cds, module_df, cell_group_df), error = function(e) {
      warnings <<- ta_append_warning(warnings, "aggregate_gene_expression_failed", paste0("Module activity by ", key, " failed: ", conditionMessage(e)))
      NULL
    })
    if (is.null(act)) {
      ta_write_tsv(ta_empty_df(c("module")), file.path(monocle_dir, file))
      return(list())
    }
    act_df <- data.frame(module = rownames(as.data.frame(act)), as.data.frame(act), row.names = NULL, check.names = FALSE)
    ta_write_tsv(act_df, file.path(monocle_dir, file))
    ta_activity_to_summary(act_df)
  }

  activity$cell_type <- write_activity(cell_type_col, "cell type", activity_files[["cell_type"]])
  activity$condition <- write_activity(condition_col, "condition", activity_files[["condition"]])
  activity$partition <- write_activity("partition", "partition", activity_files[["partition"]])
  activity$pseudotime_bin <- write_activity("pseudotime_bin", "pseudotime bin", activity_files[["pseudotime_bin"]])
  activity$branch <- write_activity("branch_path", "branch path", activity_files[["branch"]])

  modules <- lapply(split(module_df, module_df$module), function(group) {
    group$morans_I <- suppressWarnings(as.numeric(group$morans_I))
    group$q_value <- suppressWarnings(as.numeric(group$q_value))
    list(
      module = as.character(group$module[1]),
      genes_per_module = as.integer(nrow(group)),
      top_genes_by_morans_I = ta_records_from_df(group[order(group$morans_I, decreasing = TRUE), , drop = FALSE], 10),
      top_genes_by_q_value = ta_records_from_df(group[order(group$q_value), , drop = FALSE], 10),
      dominant_trend_category = if ("trend_category" %in% colnames(group) && any(!is.na(group$trend_category))) names(sort(table(group$trend_category), decreasing = TRUE))[1] else NULL
    )
  })

  list(
    gene_modules = module_df,
    gene_module_summary = list(
      num_modules = as.integer(length(modules)),
      modules = unname(modules),
      module_activity_by_cell_type = activity$cell_type,
      module_activity_by_condition = activity$condition,
      module_activity_by_partition = activity$partition,
      module_activity_by_pseudotime_bin = activity$pseudotime_bin,
      module_activity_by_branch = activity$branch
    ),
    warnings = warnings,
    confidence_flags = unique(flags)
  )
}
