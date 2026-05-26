test_that("TrajectoryAgent pseudotime QC counts finite and infinite values", {
  pseudotime <- c(cell1 = 0, cell2 = 1, cell3 = Inf, cell4 = 3, cell5 = Inf)
  partitions <- c(cell1 = "1", cell2 = "1", cell3 = "2", cell4 = "1", cell5 = "2")
  result <- ta_compute_pseudotime_qc(pseudotime, partitions, pseudotime_bins = 2)

  expect_equal(result$pseudotime_qc$finite_pseudotime_cells, 3L)
  expect_equal(result$pseudotime_qc$infinite_pseudotime_cells, 2L)
  expect_equal(result$pseudotime_qc$fraction_infinite_pseudotime, 2 / 5)
  expect_equal(result$pseudotime_qc$partitions_with_infinite_pseudotime, "2")
  expect_true(length(result$pseudotime_qc$cells_per_pseudotime_bin) > 0)
})

test_that("TrajectoryAgent pseudotime QC reports all-NA pseudotime", {
  pseudotime <- c(cell1 = NA_real_, cell2 = NA_real_)
  partitions <- c(cell1 = "1", cell2 = "1")
  result <- ta_compute_pseudotime_qc(pseudotime, partitions)

  expect_equal(result$pseudotime_qc$finite_pseudotime_cells, 0L)
  expect_equal(result$pseudotime_qc$infinite_pseudotime_cells, 0L)
  expect_true("no_pseudotime_available" %in% result$confidence_flags)
})
