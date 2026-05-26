# Trajectory-variable gene testing and directionality summaries.

ta_summarize_dynamic_gene_tables <- function(
  graph_test_df,
  dynamic_gene_trends_df = NULL,
  q_value_threshold = 0.05,
  morans_i_min = 0
) {
  if (is.null(graph_test_df) || NROW(graph_test_df) == 0) {
    graph_test_df <- ta_empty_df(c("gene_id", "gene_short_name", "q_value", "morans_I", "status"))
  }
  graph_test_df$q_value <- suppressWarnings(as.numeric(graph_test_df$q_value))
  graph_test_df$morans_I <- suppressWarnings(as.numeric(graph_test_df$morans_I))
  significant <- graph_test_df[
    !is.na(graph_test_df$q_value) &
      !is.na(graph_test_df$morans_I) &
      graph_test_df$q_value < q_value_threshold &
      graph_test_df$morans_I > morans_i_min,
    ,
    drop = FALSE
  ]
  top_moran <- significant[order(significant$morans_I, decreasing = TRUE), , drop = FALSE]
  top_q <- significant[order(significant$q_value), , drop = FALSE]

  if (is.null(dynamic_gene_trends_df) || NROW(dynamic_gene_trends_df) == 0) {
    dynamic_gene_trends_df <- ta_empty_df(c(
      "gene_id",
      "gene_short_name",
      "q_value",
      "morans_I",
      "spearman_rho",
      "spearman_p",
      "trend_category",
      "peak_pseudotime_bin"
    ))
  }
  increasing <- dynamic_gene_trends_df[as.character(dynamic_gene_trends_df$trend_category) == "increasing", , drop = FALSE]
  decreasing <- dynamic_gene_trends_df[as.character(dynamic_gene_trends_df$trend_category) == "decreasing", , drop = FALSE]
  transient <- dynamic_gene_trends_df[as.character(dynamic_gene_trends_df$trend_category) == "transient", , drop = FALSE]

  flags <- character()
  if (NROW(significant) == 0) {
    flags <- c(flags, "no_trajectory_variable_genes")
  } else if (NROW(significant) < 25) {
    flags <- c(flags, "few_trajectory_variable_genes")
  }
  if (NROW(top_moran) > 0 && is.finite(top_moran$morans_I[1]) && top_moran$morans_I[1] < 0.05) {
    flags <- c(flags, "weak_gene_dynamics")
  }

  list(
    summary = list(
      num_tested_genes = as.integer(NROW(graph_test_df)),
      num_significant_genes = as.integer(NROW(significant)),
      q_value_threshold = q_value_threshold,
      morans_i_min = morans_i_min,
      top_genes_by_morans_I = ta_records_from_df(top_moran, 20),
      top_genes_by_q_value = ta_records_from_df(top_q, 20),
      increasing_genes = ta_records_from_df(increasing[order(suppressWarnings(as.numeric(increasing$q_value))), , drop = FALSE], 20),
      decreasing_genes = ta_records_from_df(decreasing[order(suppressWarnings(as.numeric(decreasing$q_value))), , drop = FALSE], 20),
      transient_genes = ta_records_from_df(transient[order(suppressWarnings(as.numeric(transient$q_value))), , drop = FALSE], 20)
    ),
    significant_genes = significant,
    confidence_flags = unique(flags)
  )
}

ta_run_dynamic_genes <- function(
  cds,
  pseudotime,
  pseudotime_bin,
  output_dir,
  q_value_threshold = 0.05,
  morans_i_min = 0,
  top_gene_limit = 100,
  cores = 4
) {
  warnings <- list()
  graph_test_df <- tryCatch({
    result <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = cores)
    result <- as.data.frame(result, stringsAsFactors = FALSE)
    result$gene_id <- rownames(result)
    row_data <- ta_row_data(cds)
    row_data$gene_id <- rownames(row_data)
    result <- ta_normalize_gene_short_name(result, row_data, gene_id_col = "gene_id")
    keep <- intersect(c("gene_id", "gene_short_name", "q_value", "morans_I", "status"), colnames(result))
    result[, keep, drop = FALSE]
  }, error = function(e) {
    warnings <<- ta_append_warning(warnings, "graph_test_failed", paste("graph_test failed:", conditionMessage(e)))
    ta_empty_df(c("gene_id", "gene_short_name", "q_value", "morans_I", "status"))
  })
  ta_write_tsv(graph_test_df, file.path(output_dir, "monocle_outputs", "graph_test.tsv"))

  significant <- ta_summarize_dynamic_gene_tables(
    graph_test_df,
    NULL,
    q_value_threshold,
    morans_i_min
  )$significant_genes

  trend_cols <- c("gene_id", "gene_short_name", "q_value", "morans_I", "spearman_rho", "spearman_p", "trend_category", "peak_pseudotime_bin")
  dynamic_df <- ta_empty_df(trend_cols)
  finite_cells <- names(pseudotime)[is.finite(pseudotime)]
  if (NROW(significant) > 0 && length(finite_cells) >= 3) {
    significant <- significant[order(significant$q_value, -significant$morans_I), , drop = FALSE]
    trend_genes <- intersect(utils::head(significant$gene_id, top_gene_limit), rownames(cds))
    expr <- as.matrix(ta_normalized_or_counts(cds)[trend_genes, finite_cells, drop = FALSE])
    rows <- lapply(trend_genes, function(gene_id) {
      gene_expr <- as.numeric(expr[gene_id, ])
      rho <- suppressWarnings(stats::cor(gene_expr, pseudotime[finite_cells], method = "spearman", use = "complete.obs"))
      p_value <- tryCatch(
        suppressWarnings(stats::cor.test(gene_expr, pseudotime[finite_cells], method = "spearman", exact = FALSE)$p.value),
        error = function(e) NA_real_
      )
      gene_info <- significant[significant$gene_id == gene_id, , drop = FALSE][1, ]
      bin_means <- tapply(gene_expr, pseudotime_bin[finite_cells], mean, na.rm = TRUE)
      peak_bin <- if (length(bin_means) > 0) names(which.max(bin_means))[1] else NA_character_
      trend_category <- "weak_monotonic"
      if (!is.na(rho) && rho >= 0.3) {
        trend_category <- "increasing"
      } else if (!is.na(rho) && rho <= -0.3) {
        trend_category <- "decreasing"
      }
      if (length(bin_means) >= 3 && !is.na(peak_bin)) {
        peak_pos <- match(peak_bin, names(bin_means))
        if (peak_pos > 1 && peak_pos < length(bin_means) &&
          bin_means[[peak_pos]] > bin_means[[1]] &&
          bin_means[[peak_pos]] > bin_means[[length(bin_means)]]) {
          trend_category <- "transient"
        }
      }
      data.frame(
        gene_id = gene_id,
        gene_short_name = if ("gene_short_name" %in% colnames(gene_info)) gene_info$gene_short_name else gene_id,
        q_value = gene_info$q_value,
        morans_I = gene_info$morans_I,
        spearman_rho = rho,
        spearman_p = p_value,
        trend_category = trend_category,
        peak_pseudotime_bin = peak_bin,
        stringsAsFactors = FALSE
      )
    })
    dynamic_df <- do.call(rbind, rows)
  }
  ta_write_tsv(dynamic_df, file.path(output_dir, "monocle_outputs", "dynamic_gene_trends.tsv"))
  summary <- ta_summarize_dynamic_gene_tables(graph_test_df, dynamic_df, q_value_threshold, morans_i_min)
  list(
    graph_test = graph_test_df,
    dynamic_gene_trends = dynamic_df,
    significant_genes = summary$significant_genes,
    dynamic_gene_summary = summary$summary,
    warnings = warnings,
    confidence_flags = summary$confidence_flags
  )
}
