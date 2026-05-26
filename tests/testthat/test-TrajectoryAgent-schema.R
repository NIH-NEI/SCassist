test_that("TrajectoryAgent evidence JSON schema writes and LLM context excludes raw matrices", {
  evidence <- ta_create_evidence_package(
    run_config = list(seed = 1234),
    trajectory_structure = list(
      num_cells = 10L,
      num_genes = 100L,
      num_clusters = 2L,
      num_partitions = 1L,
      cells_per_partition = list("1" = 10L),
      cells_per_cluster = list("1" = 5L, "2" = 5L),
      num_branch_points = 1L,
      num_leaf_nodes = 2L,
      branch_point_ids = "Y_1",
      leaf_node_ids = c("Y_2", "Y_3")
    ),
    root_selection_evidence = list(
      method = "explicit_root_pr_nodes",
      metadata_field_used = NULL,
      metadata_value_used = NULL,
      root_pr_nodes = "Y_1",
      root_partitions = "1",
      partitions_without_roots = character(),
      candidate_root_cell_count = 1L,
      confidence = "high",
      notes = character()
    ),
    pseudotime_qc = list(
      finite_pseudotime_cells = 10L,
      infinite_pseudotime_cells = 0L,
      fraction_infinite_pseudotime = 0,
      pseudotime_range = c(0, 1),
      cells_per_pseudotime_bin = list(bin_01 = 5L, bin_02 = 5L),
      partitions_with_infinite_pseudotime = character()
    ),
    branch_summary = list(),
    dynamic_gene_summary = ta_empty_dynamic_gene_summary(),
    gene_module_summary = ta_empty_gene_module_summary(),
    metadata_associations = ta_empty_metadata_associations(),
    confidence_flags = character(),
    confidence_score = 0.9,
    confidence_label = "high"
  )
  evidence$llm_interpretation_input <- ta_deterministic_interpretation(evidence)
  ta_validate_evidence_schema(evidence)

  out <- tempfile(fileext = ".json")
  jsonlite::write_json(evidence, out, auto_unbox = TRUE, pretty = TRUE, na = "null")
  expect_true(file.exists(out))
  parsed <- jsonlite::read_json(out)
  expect_true(all(c("trajectory_structure", "pseudotime_qc", "confidence_score") %in% names(parsed)))

  llm_context <- ta_build_llm_context(evidence)
  ta_validate_llm_context(llm_context)
  expect_false("expression_matrix" %in% names(llm_context))
  expect_false("graph_test" %in% names(llm_context))
  expect_equal(
    sort(names(llm_context)),
    sort(c(
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
      "cautious_interpretation_rules"
    ))
  )
})

test_that("TrajectoryAgent LLM context validator rejects raw and oversized content", {
  expect_error(
    ta_validate_llm_context(list(expression_matrix = matrix(1, nrow = 2, ncol = 2))),
    "forbidden raw field"
  )
  expect_error(
    ta_validate_llm_context(list(summary = data.frame(x = seq_len(101))), max_items_per_section = 100),
    "oversized data frame"
  )
  compact <- list(dynamic_gene_summary = list(top_genes_by_q_value = list(list(gene_id = "g1"))))
  expect_silent(ta_validate_llm_context(compact))
})

test_that("TrajectoryAgent report avoids prohibited overclaiming phrases", {
  evidence <- ta_create_evidence_package(
    run_config = list(seed = 1234),
    trajectory_structure = list(
      num_cells = 10L,
      num_genes = 100L,
      num_clusters = 2L,
      num_partitions = 1L,
      cells_per_partition = list("1" = 10L),
      cells_per_cluster = list("1" = 5L),
      num_branch_points = 0L,
      num_leaf_nodes = 0L,
      branch_point_ids = character(),
      leaf_node_ids = character()
    ),
    root_selection_evidence = list(
      method = "explicit_root_pr_nodes",
      metadata_field_used = NULL,
      metadata_value_used = NULL,
      root_pr_nodes = "Y_1",
      root_partitions = "1",
      partitions_without_roots = character(),
      candidate_root_cell_count = 1L,
      confidence = "high",
      notes = character()
    ),
    pseudotime_qc = list(
      finite_pseudotime_cells = 10L,
      infinite_pseudotime_cells = 0L,
      fraction_infinite_pseudotime = 0,
      pseudotime_range = c(0, 1),
      cells_per_pseudotime_bin = list(bin_01 = 5L, bin_02 = 5L),
      partitions_with_infinite_pseudotime = character()
    ),
    branch_summary = list(),
    dynamic_gene_summary = within(ta_empty_dynamic_gene_summary(), {
      num_tested_genes <- 3L
      num_significant_genes <- 1L
      top_genes_by_q_value <- list(list(gene_id = "g1", gene_short_name = "G1", q_value = 0.01, morans_I = 0.2))
      top_genes_by_morans_I <- top_genes_by_q_value
    }),
    gene_module_summary = ta_empty_gene_module_summary(),
    metadata_associations = ta_empty_metadata_associations(),
    confidence_flags = character(),
    confidence_score = 0.8,
    confidence_label = "high"
  )
  evidence$llm_interpretation_input <- ta_deterministic_interpretation(evidence)
  interpretation_text <- tolower(paste(unlist(evidence$llm_interpretation_input), collapse = "\n"))
  output_dir <- tempfile("trajectory-report-")
  ta_write_trajectory_reports(evidence, output_dir)
  report <- paste(readLines(file.path(output_dir, "trajectory_report.md")), collapse = "\n")
  report_lower <- tolower(report)
  prohibited <- c("proves", "causes", "definitively represents", "confirms lineage commitment", "Monocle proves", "gene X causes")
  expect_false(any(vapply(tolower(prohibited), grepl, logical(1), x = report_lower, fixed = TRUE)))
  expect_false(any(vapply(tolower(prohibited), grepl, logical(1), x = interpretation_text, fixed = TRUE)))
  expect_true(grepl("consistent with", report_lower, fixed = TRUE))
  expect_true(grepl("suggesting", report_lower, fixed = TRUE))
  expect_true(grepl("may reflect", report_lower, fixed = TRUE))
  expect_true(grepl("supports, but does not prove", report_lower, fixed = TRUE))
})
