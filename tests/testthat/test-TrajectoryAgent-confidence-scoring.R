test_that("TrajectoryAgent confidence scoring applies deterministic penalties", {
  root <- list(
    root_pr_nodes = character(),
    partitions_without_roots = "1"
  )
  qc <- list(fraction_infinite_pseudotime = 0.3)
  dynamic <- list(num_significant_genes = 0L)
  modules <- list(num_modules = 0L)
  metadata <- list(possible_batch_or_donor_confounding = "batch_strongly_associated_with_pseudotime")
  branches <- list()

  result <- ta_score_trajectory_confidence(root, qc, dynamic, modules, metadata, branches)
  expect_equal(result$label, "low")
  expect_true(result$score < 0.45)
  expect_true("no_defensible_root_selected" %in% result$flags)
  expect_true("high_infinite_pseudotime_fraction" %in% result$flags)
  expect_true("no_trajectory_variable_genes" %in% result$flags)
  expect_true("batch_or_donor_strongly_associated_with_pseudotime" %in% result$flags)
})
