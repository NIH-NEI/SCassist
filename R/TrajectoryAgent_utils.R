# Utility helpers for the R-native SCAssist TrajectoryAgent.

ta_null_arg <- function(x) {
  is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

ta_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "TrajectoryAgent requires missing R package(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

ta_warning <- function(flag, message, details = list()) {
  list(flag = flag, message = message, details = details)
}

ta_append_warning <- function(warnings, flag, message, details = list()) {
  c(warnings, list(ta_warning(flag, message, details)))
}

ta_as_list <- function(x) {
  if (is.null(x)) {
    return(list())
  }
  if (is.list(x) && is.null(dim(x))) {
    return(x)
  }
  as.list(x)
}

ta_named_count_list <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  counts <- table(x)
  out <- as.list(as.integer(counts))
  names(out) <- names(counts)
  out
}

ta_records_from_df <- function(df, n = 20) {
  if (is.null(df) || NROW(df) == 0) {
    return(list())
  }
  df <- as.data.frame(utils::head(df, n), stringsAsFactors = FALSE)
  lapply(seq_len(nrow(df)), function(i) as.list(df[i, , drop = FALSE]))
}

ta_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    x,
    path,
    pretty = TRUE,
    auto_unbox = TRUE,
    na = "null",
    null = "null"
  )
}

ta_write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    as.data.frame(df, stringsAsFactors = FALSE),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

ta_empty_df <- function(columns) {
  as.data.frame(setNames(rep(list(character()), length(columns)), columns))
}

ta_first_nonempty <- function(...) {
  values <- list(...)
  if (length(values) == 0) {
    return(NULL)
  }
  n <- max(vapply(values, length, integer(1)))
  out <- rep(NA_character_, n)
  for (value in values) {
    if (is.null(value)) {
      next
    }
    value <- as.character(value)
    if (length(value) == 1 && n > 1) {
      value <- rep(value, n)
    }
    value <- value[seq_len(n)]
    take <- (is.na(out) | !nzchar(out)) & !is.na(value) & nzchar(value)
    out[take] <- value[take]
  }
  out
}

ta_normalize_gene_short_name <- function(df, gene_metadata = NULL, gene_id_col = "gene_id") {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!gene_id_col %in% colnames(df)) {
    stop("gene_id column is required for gene_short_name normalization.", call. = FALSE)
  }
  original_order <- seq_len(nrow(df))
  df$.ta_original_order <- original_order

  gene_metadata_col <- NULL
  if (!is.null(gene_metadata)) {
    gene_metadata <- as.data.frame(gene_metadata, stringsAsFactors = FALSE)
    if (!gene_id_col %in% colnames(gene_metadata)) {
      gene_metadata[[gene_id_col]] <- rownames(gene_metadata)
    }
    if ("gene_short_name" %in% colnames(gene_metadata)) {
      gene_metadata <- gene_metadata[, c(gene_id_col, "gene_short_name"), drop = FALSE]
      colnames(gene_metadata)[colnames(gene_metadata) == "gene_short_name"] <- ".ta_gene_short_name_metadata"
      df <- merge(df, gene_metadata, by = gene_id_col, all.x = TRUE, sort = FALSE)
    }
  }
  if (".ta_gene_short_name_metadata" %in% colnames(df)) {
    gene_metadata_col <- df$.ta_gene_short_name_metadata
  }

  existing <- if ("gene_short_name" %in% colnames(df)) df$gene_short_name else NULL
  existing_x <- if ("gene_short_name.x" %in% colnames(df)) df$gene_short_name.x else NULL
  existing_y <- if ("gene_short_name.y" %in% colnames(df)) df$gene_short_name.y else NULL
  fallback <- df[[gene_id_col]]
  df$gene_short_name <- ta_first_nonempty(existing, existing_x, existing_y, gene_metadata_col, fallback)

  remove_cols <- intersect(c("gene_short_name.x", "gene_short_name.y", ".ta_gene_short_name_metadata"), colnames(df))
  if (length(remove_cols) > 0) {
    df <- df[, setdiff(colnames(df), remove_cols), drop = FALSE]
  }
  df <- df[order(df$.ta_original_order), , drop = FALSE]
  df$.ta_original_order <- NULL
  rownames(df) <- NULL
  df
}

ta_check_requested_metadata_columns <- function(colData_df, requested_cols) {
  colData_df <- as.data.frame(colData_df, stringsAsFactors = FALSE)
  requested_cols <- unique(requested_cols[!vapply(requested_cols, ta_null_arg, logical(1))])
  present <- requested_cols[requested_cols %in% colnames(colData_df)]
  missing <- requested_cols[!requested_cols %in% colnames(colData_df)]
  warnings <- lapply(missing, function(column) {
    ta_warning(
      paste0("requested_metadata_column_missing:", column),
      paste0("Requested metadata column '", column, "' is absent and will be skipped."),
      list(column = column)
    )
  })
  list(present = present, missing = missing, warnings = warnings)
}

ta_formula_missing_columns <- function(formula_string, colData_df) {
  if (ta_null_arg(formula_string)) {
    return(character())
  }
  vars <- tryCatch(all.vars(stats::as.formula(paste("~", formula_string))), error = function(e) character())
  setdiff(vars, colnames(colData_df))
}

ta_parse_root_cells_arg <- function(root_cells) {
  if (ta_null_arg(root_cells)) {
    return(NULL)
  }
  if (length(root_cells) == 1 && file.exists(root_cells)) {
    table <- tryCatch(
      utils::read.delim(root_cells, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (!is.null(table) && "cell_id" %in% colnames(table)) {
      return(as.character(table$cell_id[!is.na(table$cell_id) & nzchar(as.character(table$cell_id))]))
    }
    cells <- readLines(root_cells, warn = FALSE)
    return(cells[!is.na(cells) & nzchar(cells)])
  }
  trimws(unlist(strsplit(as.character(root_cells), ",", fixed = TRUE), use.names = FALSE))
}

ta_validate_llm_context <- function(llm_context, max_items_per_section = 100, max_json_chars = 200000) {
  forbidden_names <- c(
    "expression_matrix",
    "counts",
    "matrix",
    "raw_counts",
    "normalized_counts",
    "X",
    "assays",
    "cell_trajectory",
    "umap",
    "graph_test",
    "module_activity_matrix"
  )
  json <- jsonlite::toJSON(llm_context, auto_unbox = TRUE, na = "null", null = "null")
  if (nchar(json, type = "chars") > max_json_chars) {
    stop("llm_context is too large to be LLM-ready.", call. = FALSE)
  }
  check_node <- function(x, path = character()) {
    if (length(path) > 0 && tail(path, 1) %in% forbidden_names) {
      stop("llm_context contains forbidden raw field: ", paste(path, collapse = "."), call. = FALSE)
    }
    if (is.data.frame(x)) {
      if (nrow(x) > max_items_per_section) {
        stop("llm_context contains oversized data frame: ", paste(path, collapse = "."), call. = FALSE)
      }
      return(invisible(TRUE))
    }
    if (is.matrix(x) || inherits(x, "Matrix")) {
      stop("llm_context contains matrix-like raw data: ", paste(path, collapse = "."), call. = FALSE)
    }
    if (is.list(x)) {
      if (length(x) > max_items_per_section && is.null(names(x))) {
        stop("llm_context contains oversized unnamed list: ", paste(path, collapse = "."), call. = FALSE)
      }
      nms <- names(x)
      for (i in seq_along(x)) {
        child_name <- if (!is.null(nms) && nzchar(nms[[i]])) nms[[i]] else paste0("[", i, "]")
        check_node(x[[i]], c(path, child_name))
      }
    }
    invisible(TRUE)
  }
  check_node(llm_context)
  invisible(TRUE)
}

ta_composition_summary <- function(values, top_n = 5) {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) {
    return(list())
  }
  counts <- sort(table(values), decreasing = TRUE)
  total <- sum(counts)
  lapply(seq_len(min(length(counts), top_n)), function(i) {
    list(
      value = names(counts)[i],
      count = as.integer(counts[[i]]),
      fraction = as.numeric(counts[[i]] / total)
    )
  })
}

ta_make_pseudotime_bins <- function(pseudotime, n_bins = 10) {
  bins <- rep(NA_character_, length(pseudotime))
  names(bins) <- names(pseudotime)
  finite_idx <- which(is.finite(pseudotime))
  if (length(finite_idx) == 0) {
    return(bins)
  }
  finite_pt <- pseudotime[finite_idx]
  if (length(unique(finite_pt)) < 2) {
    bins[finite_idx] <- "bin_01"
    return(bins)
  }
  n_bins <- max(1, min(n_bins, length(unique(finite_pt))))
  breaks <- unique(as.numeric(stats::quantile(
    finite_pt,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE,
    type = 7
  )))
  if (length(breaks) < 2) {
    bins[finite_idx] <- "bin_01"
    return(bins)
  }
  labels <- sprintf("bin_%02d", seq_len(length(breaks) - 1))
  bins[finite_idx] <- as.character(cut(
    finite_pt,
    breaks = breaks,
    include.lowest = TRUE,
    labels = labels
  ))
  bins
}

ta_col_data <- function(cds) {
  if (is.list(cds) && !is.null(cds$cell_metadata)) {
    return(as.data.frame(cds$cell_metadata, stringsAsFactors = FALSE))
  }
  as.data.frame(SummarizedExperiment::colData(cds), stringsAsFactors = FALSE)
}

ta_row_data <- function(cds) {
  if (is.list(cds) && !is.null(cds$gene_metadata)) {
    return(as.data.frame(cds$gene_metadata, stringsAsFactors = FALSE))
  }
  as.data.frame(SummarizedExperiment::rowData(cds), stringsAsFactors = FALSE)
}

ta_safe_col <- function(df, column) {
  if (ta_null_arg(column) || !column %in% colnames(df)) {
    return(NULL)
  }
  df[[column]]
}

ta_normalized_or_counts <- function(cds) {
  norm <- tryCatch(monocle3::normalized_counts(cds), error = function(e) NULL)
  if (!is.null(norm)) {
    return(norm)
  }
  SummarizedExperiment::assay(cds)
}

ta_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(NULL)
  }
  as.character(utils::packageVersion(pkg))
}
