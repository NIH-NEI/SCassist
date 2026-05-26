test_that("TrajectoryAgent skips gene modules when significant genes are below threshold", {
  output_dir <- tempfile("trajectory-modules-")
  result <- ta_run_gene_modules(
    cds = NULL,
    significant_genes = data.frame(gene_id = paste0("g", 1:9)),
    graph_test_df = data.frame(gene_id = paste0("g", 1:9), q_value = 0.01, morans_I = 0.1),
    dynamic_gene_trends_df = data.frame(),
    cell_trajectory = data.frame(cell_id = character()),
    output_dir = output_dir
  )

  expect_true("too_few_genes_for_gene_modules" %in% result$confidence_flags)
  expect_true(file.exists(file.path(output_dir, "monocle_outputs", "gene_modules.tsv")))
  expect_equal(nrow(utils::read.delim(file.path(output_dir, "monocle_outputs", "gene_modules.tsv"))), 0)
})

test_that("TrajectoryAgent module helper attempts find_gene_modules when threshold is met", {
  module_file <- file.path("R", "TrajectoryAgent_gene_modules.R")
  if (!file.exists(module_file)) {
    module_file <- file.path("..", "..", "R", "TrajectoryAgent_gene_modules.R")
  }
  files <- readLines(module_file, warn = FALSE)
  text <- paste(files, collapse = "\n")
  expect_true(grepl("NROW\\(significant_genes\\) < 10", text))
  expect_true(grepl("find_gene_modules\\(", text))
  expect_true(grepl("aggregate_gene_expression\\(", text))
})
