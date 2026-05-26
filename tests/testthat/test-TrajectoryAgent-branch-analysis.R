test_that("TrajectoryAgent detects branch and leaf nodes and flags small paths", {
  skip_if_not_installed("igraph")
  graph <- igraph::make_graph(c("Y_1", "Y_2", "Y_1", "Y_3", "Y_1", "Y_4"), directed = FALSE)
  closest <- c(
    cell1 = "Y_2",
    cell2 = "Y_2",
    cell3 = "Y_3",
    cell4 = "Y_4",
    cell5 = "Y_1"
  )
  cell_trajectory <- data.frame(
    cell_id = names(closest),
    pseudotime = c(1, 2, 3, 4, 0),
    pseudotime_bin = c("bin_01", "bin_01", "bin_02", "bin_02", "bin_01"),
    closest_pr_node = unname(closest),
    cell_type = c("A", "A", "B", "C", "A"),
    stringsAsFactors = FALSE
  )
  result <- ta_analyze_branches(
    graph,
    closest,
    cell_trajectory,
    cell_type_col = "cell_type",
    min_cells_per_branch_path = 3
  )

  expect_equal(length(result$branch_summary), 1)
  expect_equal(result$branch_summary[[1]]$branch_node_id, "Y_1")
  expect_true("branch_paths_below_min_cell_count" %in% result$confidence_flags)
  expect_equal(sort(igraph::V(graph)$name[igraph::degree(graph) == 1]), c("Y_2", "Y_3", "Y_4"))
})

test_that("TrajectoryAgent homogeneous adequate branch paths are low confidence", {
  skip_if_not_installed("igraph")
  graph <- igraph::make_graph(c("Y_1", "Y_2", "Y_1", "Y_3", "Y_1", "Y_4"), directed = FALSE)
  closest <- c(
    setNames(rep("Y_2", 4), paste0("a", 1:4)),
    setNames(rep("Y_3", 4), paste0("b", 1:4)),
    setNames(rep("Y_4", 4), paste0("c", 1:4))
  )
  cell_trajectory <- data.frame(
    cell_id = names(closest),
    pseudotime = seq_along(closest),
    pseudotime_bin = "bin_01",
    closest_pr_node = unname(closest),
    cell_type = "same",
    condition = "same_condition",
    stringsAsFactors = FALSE
  )
  result <- ta_analyze_branches(
    graph,
    closest,
    cell_trajectory,
    cell_type_col = "cell_type",
    condition_col = "condition",
    min_cells_per_branch_path = 3
  )
  branch <- result$branch_summary[[1]]
  expect_equal(branch$confidence, "low")
  expect_true("branch_paths_not_distinguishable" %in% branch$confidence_flags)
  expect_true("homogeneous_branch_paths" %in% branch$confidence_flags)
})

test_that("TrajectoryAgent distinguishable branch paths without gene support are not high", {
  skip_if_not_installed("igraph")
  graph <- igraph::make_graph(c("Y_1", "Y_2", "Y_1", "Y_3", "Y_1", "Y_4"), directed = FALSE)
  closest <- c(
    setNames(rep("Y_2", 4), paste0("a", 1:4)),
    setNames(rep("Y_3", 4), paste0("b", 1:4)),
    setNames(rep("Y_4", 4), paste0("c", 1:4))
  )
  cell_trajectory <- data.frame(
    cell_id = names(closest),
    pseudotime = seq_along(closest),
    pseudotime_bin = "bin_01",
    closest_pr_node = unname(closest),
    cell_type = rep(c("A", "B", "C"), each = 4),
    stringsAsFactors = FALSE
  )
  result <- ta_analyze_branches(
    graph,
    closest,
    cell_trajectory,
    cell_type_col = "cell_type",
    min_cells_per_branch_path = 3
  )
  branch <- result$branch_summary[[1]]
  expect_true(branch$confidence %in% c("medium", "low"))
  expect_false(identical(branch$confidence, "high"))
  expect_true("no_branch_associated_genes_or_modules" %in% branch$confidence_flags)
})

test_that("TrajectoryAgent supported distinguishable branches can be high confidence", {
  skip_if_not_installed("igraph")
  graph <- igraph::make_graph(c("Y_1", "Y_2", "Y_1", "Y_3", "Y_1", "Y_4"), directed = FALSE)
  closest <- c(
    setNames(rep("Y_2", 4), paste0("a", 1:4)),
    setNames(rep("Y_3", 4), paste0("b", 1:4)),
    setNames(rep("Y_4", 4), paste0("c", 1:4))
  )
  cell_trajectory <- data.frame(
    cell_id = names(closest),
    pseudotime = seq_along(closest),
    pseudotime_bin = "bin_01",
    closest_pr_node = unname(closest),
    cell_type = rep(c("A", "B", "C"), each = 4),
    stringsAsFactors = FALSE
  )
  result <- ta_analyze_branches(
    graph,
    closest,
    cell_trajectory,
    cell_type_col = "cell_type",
    min_cells_per_branch_path = 3,
    branch_associated_genes = list(Y_1 = list(list(gene_id = "g1")))
  )
  expect_equal(result$branch_summary[[1]]$confidence, "high")
})

test_that("TrajectoryAgent branches with one adequate path are low confidence", {
  skip_if_not_installed("igraph")
  graph <- igraph::make_graph(c("Y_1", "Y_2", "Y_1", "Y_3", "Y_1", "Y_4"), directed = FALSE)
  closest <- c(setNames(rep("Y_2", 5), paste0("a", 1:5)), b1 = "Y_3", c1 = "Y_4")
  cell_trajectory <- data.frame(
    cell_id = names(closest),
    pseudotime = seq_along(closest),
    pseudotime_bin = "bin_01",
    closest_pr_node = unname(closest),
    cell_type = c(rep("A", 5), "B", "C"),
    stringsAsFactors = FALSE
  )
  result <- ta_analyze_branches(
    graph,
    closest,
    cell_trajectory,
    cell_type_col = "cell_type",
    min_cells_per_branch_path = 3
  )
  expect_equal(result$branch_summary[[1]]$confidence, "low")
  expect_true("no_adequate_branch_paths" %in% result$branch_summary[[1]]$confidence_flags)
})
