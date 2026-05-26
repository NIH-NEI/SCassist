test_that("TrajectoryAgent does not contain interactive Monocle3 calls or Python references", {
  r_dir <- if (dir.exists("R")) "R" else file.path("..", "..", "R")
  files <- list.files(r_dir, pattern = "^TrajectoryAgent_.*\\.R$", full.names = TRUE)
  text <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
  interactive_calls <- paste0(c("choose", "choose"), c("_cells", "_graph_segments"))
  for (call_name in interactive_calls) {
    expect_false(grepl(paste0(call_name, "\\s*\\("), text))
  }
  forbidden_python_terms <- paste(
    c(
      paste0("reti", "culate"),
      paste0("ann", "data"),
      paste0("scan", "py"),
      paste0("pan", "das"),
      paste0("num", "py"),
      paste0("sci", "py"),
      paste0("pyd", "antic")
    ),
    collapse = "|"
  )
  expect_false(grepl(forbidden_python_terms, text, ignore.case = TRUE))
})

test_that("TrajectoryAgent order_cells calls are guarded by explicit roots", {
  r_dir <- if (dir.exists("R")) "R" else file.path("..", "..", "R")
  files <- list.files(r_dir, pattern = "^TrajectoryAgent_.*\\.R$", full.names = TRUE)
  lines <- unlist(lapply(files, readLines, warn = FALSE))
  order_lines <- lines[grepl("order_cells\\s*\\(", lines)]
  expect_true(length(order_lines) > 0)
  expect_true(all(grepl("root_pr_nodes\\s*=", order_lines) | grepl("root_cells\\s*=", order_lines)))
  expect_false(any(grepl("order_cells\\s*\\(\\s*cds\\s*\\)", order_lines)))
})

test_that("TrajectoryAgent graph_test uses principal_graph", {
  r_dir <- if (dir.exists("R")) "R" else file.path("..", "..", "R")
  files <- list.files(r_dir, pattern = "^TrajectoryAgent_.*\\.R$", full.names = TRUE)
  text <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
  expect_true(grepl("graph_test\\s*\\(", text))
  expect_true(grepl("neighbor_graph\\s*=\\s*\"principal_graph\"", text))
})
