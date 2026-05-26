test_that("TrajectoryAgent optional Monocle3 integration produces core reports", {
  skip_if_not_installed("monocle3")
  skip_if_not_installed("Matrix")

  set.seed(1)
  counts <- Matrix::Matrix(matrix(rpois(30 * 20, lambda = 3), nrow = 30), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  cell_metadata <- data.frame(
    timepoint = rep(c("day0", "day1"), each = 10),
    cell_type = rep(c("start", "late"), each = 10),
    row.names = colnames(counts)
  )
  gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  output_dir <- tempfile("trajectory-agent-")

  result <- run_trajectory_agent(
    expression_matrix = counts,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata,
    output_dir = output_dir,
    time_col = "timepoint",
    earliest_time_value = "day0",
    num_dim = 5,
    min_cells_per_partition = 5,
    min_root_candidate_cells = 2,
    cores = 1,
    save_cds = FALSE
  )

  expect_equal(result$agent_name, "TrajectoryAgent")
  expect_true(file.exists(file.path(output_dir, "trajectory_evidence.json")))
  expect_true(file.exists(file.path(output_dir, "llm_context.json")))
  expect_true(file.exists(file.path(output_dir, "trajectory_report.md")))
})
