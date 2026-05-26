test_that("TrajectoryAgent dynamic gene summary sorts top genes and trend categories", {
  graph_test <- data.frame(
    gene_id = c("g1", "g2", "g3"),
    gene_short_name = c("G1", "G2", "G3"),
    q_value = c(0.001, 0.02, 0.2),
    morans_I = c(0.4, 0.2, 0.9),
    stringsAsFactors = FALSE
  )
  trends <- data.frame(
    gene_id = c("g1", "g2", "g4"),
    gene_short_name = c("G1", "G2", "G4"),
    q_value = c(0.001, 0.02, 0.01),
    morans_I = c(0.4, 0.2, 0.3),
    trend_category = c("increasing", "decreasing", "transient"),
    stringsAsFactors = FALSE
  )
  result <- ta_summarize_dynamic_gene_tables(graph_test, trends, q_value_threshold = 0.05)

  expect_equal(result$summary$num_significant_genes, 2L)
  expect_equal(result$summary$top_genes_by_morans_I[[1]]$gene_id, "g1")
  expect_equal(result$summary$top_genes_by_q_value[[1]]$gene_id, "g1")
  expect_equal(result$summary$increasing_genes[[1]]$gene_id, "g1")
  expect_equal(result$summary$decreasing_genes[[1]]$gene_id, "g2")
  expect_equal(result$summary$transient_genes[[1]]$gene_id, "g4")
})

test_that("TrajectoryAgent normalizes graph_test gene_short_name columns", {
  graph_test <- data.frame(
    gene_id = c("g1", "g2"),
    gene_short_name = c("Existing1", ""),
    stringsAsFactors = FALSE
  )
  gene_metadata <- data.frame(
    gene_id = c("g1", "g2"),
    gene_short_name = c("Meta1", "Meta2"),
    stringsAsFactors = FALSE
  )
  result <- ta_normalize_gene_short_name(graph_test, gene_metadata)
  expect_equal(sum(colnames(result) == "gene_short_name"), 1)
  expect_false(any(c("gene_short_name.x", "gene_short_name.y") %in% colnames(result)))
  expect_equal(result$gene_short_name, c("Existing1", "Meta2"))
  expect_equal(result$gene_id, c("g1", "g2"))
})

test_that("TrajectoryAgent gene_short_name normalization coalesces suffixes and falls back to gene_id", {
  module_table <- data.frame(
    gene_id = c("g1", "g2", "g3"),
    gene_short_name.x = c("", "X2", ""),
    gene_short_name.y = c("Y1", "", ""),
    stringsAsFactors = FALSE
  )
  result <- ta_normalize_gene_short_name(module_table)
  expect_equal(sum(colnames(result) == "gene_short_name"), 1)
  expect_false(any(c("gene_short_name.x", "gene_short_name.y") %in% colnames(result)))
  expect_equal(result$gene_short_name, c("Y1", "X2", "g3"))
})
