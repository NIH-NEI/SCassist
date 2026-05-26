test_that("TrajectoryAgent CLI root_cells helper parses comma-separated values", {
  expect_equal(ta_parse_root_cells_arg("cell1, cell2,cell3"), c("cell1", "cell2", "cell3"))
})

test_that("TrajectoryAgent CLI root_cells helper parses one-cell-per-line files", {
  path <- tempfile()
  writeLines(c("cell1", "cell2"), path)
  expect_equal(ta_parse_root_cells_arg(path), c("cell1", "cell2"))
})

test_that("TrajectoryAgent CLI root_cells helper parses cell_id TSV files", {
  path <- tempfile(fileext = ".tsv")
  utils::write.table(
    data.frame(cell_id = c("cellA", "cellB"), extra = 1:2),
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  expect_equal(ta_parse_root_cells_arg(path), c("cellA", "cellB"))
})
