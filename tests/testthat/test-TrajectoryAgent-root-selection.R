test_that("TrajectoryAgent explicit root principal graph nodes are honored", {
  skip_if_not_installed("igraph")
  graph <- igraph::make_graph(c("Y_1", "Y_2", "Y_2", "Y_3"), directed = FALSE)
  mock_cds <- list(
    graph = graph,
    closest_pr_vertex = c(cell1 = "Y_1", cell2 = "Y_2", cell3 = "Y_3"),
    partitions = c(cell1 = "1", cell2 = "1", cell3 = "1"),
    cell_metadata = data.frame(row.names = c("cell1", "cell2", "cell3"))
  )
  result <- ta_select_roots(
    mock_cds,
    root_pr_nodes = c("Y_2", "missing"),
    min_cells_per_partition = 1,
    min_root_candidate_cells = 1
  )
  expect_equal(result$root_pr_nodes, "Y_2")
  expect_equal(result$method, "explicit_root_pr_nodes")
})

test_that("TrajectoryAgent candidate cells choose most frequent graph node per partition", {
  skip_if_not_installed("igraph")
  graph <- igraph::make_graph(c("Y_1", "Y_2", "Y_2", "Y_3", "Y_4", "Y_5"), directed = FALSE)
  mock_cds <- list(
    graph = graph,
    closest_pr_vertex = c(cell1 = "Y_1", cell2 = "Y_1", cell3 = "Y_2", cell4 = "Y_4", cell5 = "Y_4", cell6 = "Y_5"),
    partitions = c(cell1 = "1", cell2 = "1", cell3 = "1", cell4 = "2", cell5 = "2", cell6 = "2")
  )
  result <- ta_select_roots_from_candidate_cells(
    mock_cds,
    candidate_cell_ids = c("cell1", "cell2", "cell4", "cell5"),
    min_cells_per_partition = 1,
    min_root_candidate_cells = 2
  )
  expect_equal(sort(result$root_pr_nodes), c("Y_1", "Y_4"))
  expect_equal(sort(result$root_partitions), c("1", "2"))
})

test_that("TrajectoryAgent partitions without candidate roots are flagged", {
  skip_if_not_installed("igraph")
  graph <- igraph::make_graph(c("Y_1", "Y_2", "Y_3", "Y_4"), directed = FALSE)
  mock_cds <- list(
    graph = graph,
    closest_pr_vertex = c(cell1 = "Y_1", cell2 = "Y_2", cell3 = "Y_3", cell4 = "Y_4"),
    partitions = c(cell1 = "1", cell2 = "1", cell3 = "2", cell4 = "2")
  )
  result <- ta_select_roots_from_candidate_cells(
    mock_cds,
    candidate_cell_ids = c("cell1", "cell2"),
    min_cells_per_partition = 1,
    min_root_candidate_cells = 2
  )
  expect_equal(result$root_pr_nodes, "Y_1")
  expect_equal(result$partitions_without_roots, "2")
})
