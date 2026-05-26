#!/usr/bin/env Rscript

parse_args <- function(argv) {
  args <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    name <- sub("^--", "", key)
    value <- TRUE
    if (i + 1 <= length(argv) && !startsWith(argv[[i + 1]], "--")) {
      value <- argv[[i + 1]]
      i <- i + 1
    }
    args[[name]] <- value
    i <- i + 1
  }
  args
}

arg_or_null <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || identical(value, "")) NULL else value
}

as_num <- function(args, name, default) {
  value <- arg_or_null(args, name)
  if (is.null(value)) default else as.numeric(value)
}

as_int <- function(args, name, default) {
  as.integer(as_num(args, name, default))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$output_dir)) {
  stop("--output_dir is required", call. = FALSE)
}

loaded_from_package <- FALSE
if (requireNamespace("SCassist", quietly = TRUE)) {
  suppressPackageStartupMessages(library(SCassist))
  loaded_from_package <- exists("run_trajectory_agent", mode = "function")
}
if (!loaded_from_package) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else "inst/scripts/run_trajectory_agent.R"
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
  r_files <- list.files(file.path(repo_root, "R"), pattern = "^TrajectoryAgent_.*\\.R$", full.names = TRUE)
  for (file in sort(r_files)) source(file)
}

cds <- if (!is.null(args$cds_rds)) readRDS(args$cds_rds) else NULL
matrix <- if (!is.null(args$matrix_rds)) readRDS(args$matrix_rds) else NULL
cell_metadata <- if (!is.null(args$cell_metadata_tsv)) {
  df <- utils::read.delim(args$cell_metadata_tsv, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- df[[1]]
  df[, -1, drop = FALSE]
} else NULL
gene_metadata <- if (!is.null(args$gene_metadata_tsv)) {
  df <- utils::read.delim(args$gene_metadata_tsv, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- df[[1]]
  df[, -1, drop = FALSE]
} else NULL
seurat_object <- if (!is.null(args$seurat_rds)) readRDS(args$seurat_rds) else NULL
root_pr_nodes <- arg_or_null(args, "root_pr_nodes")
if (!is.null(root_pr_nodes)) {
  root_pr_nodes <- trimws(strsplit(root_pr_nodes, ",", fixed = TRUE)[[1]])
}
root_cells <- ta_parse_root_cells_arg(arg_or_null(args, "root_cells"))

run_trajectory_agent(
  cds = cds,
  expression_matrix = matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata,
  seurat_object = seurat_object,
  output_dir = args$output_dir,
  num_dim = as_int(args, "num_dim", 50),
  alignment_group = arg_or_null(args, "alignment_group"),
  residual_model_formula_str = arg_or_null(args, "residual_model_formula_str"),
  cell_type_col = arg_or_null(args, "cell_type_col"),
  condition_col = arg_or_null(args, "condition_col"),
  batch_col = arg_or_null(args, "batch_col"),
  donor_col = arg_or_null(args, "donor_col"),
  time_col = arg_or_null(args, "time_col"),
  earliest_time_value = arg_or_null(args, "earliest_time_value"),
  root_cell_type_col = arg_or_null(args, "root_cell_type_col"),
  root_cell_type_value = arg_or_null(args, "root_cell_type_value"),
  root_condition_col = arg_or_null(args, "root_condition_col"),
  root_condition_value = arg_or_null(args, "root_condition_value"),
  root_pr_nodes = root_pr_nodes,
  root_cells = root_cells,
  q_value_threshold = as_num(args, "q_value_threshold", 0.05),
  morans_i_min = as_num(args, "morans_i_min", 0),
  top_gene_limit = as_int(args, "top_gene_limit", 100),
  min_cells_per_partition = as_int(args, "min_cells_per_partition", 30),
  min_cells_per_branch_path = as_int(args, "min_cells_per_branch_path", 20),
  pseudotime_bins = as_int(args, "pseudotime_bins", 10),
  cores = as_int(args, "cores", 4),
  seed = as_int(args, "seed", 1234)
)
