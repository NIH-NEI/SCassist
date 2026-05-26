test_that("TrajectoryAgent missing optional metadata columns produce warnings", {
  col_data <- data.frame(existing = c("a", "b"), row.names = c("cell1", "cell2"))
  result <- ta_check_requested_metadata_columns(col_data, c("cell_type", "existing"))

  expect_equal(result$present, "existing")
  expect_equal(result$missing, "cell_type")
  expect_equal(result$warnings[[1]]$flag, "requested_metadata_column_missing:cell_type")
})

test_that("TrajectoryAgent missing alignment group is represented as a warning flag", {
  col_data <- data.frame(batch = c("b1", "b2"))
  result <- ta_check_requested_metadata_columns(col_data, c("missing_alignment"))

  expect_equal(result$warnings[[1]]$flag, "requested_metadata_column_missing:missing_alignment")
})

test_that("TrajectoryAgent missing root metadata column can be flagged without crashing", {
  col_data <- data.frame(cell_type = c("A", "B"))
  result <- ta_check_requested_metadata_columns(col_data, c("missing_root_col"))

  expect_equal(result$missing, "missing_root_col")
  expect_equal(result$warnings[[1]]$flag, "requested_metadata_column_missing:missing_root_col")
})
