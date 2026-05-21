#' @title Compare CellChat Communication Between Two Conditions
#'
#' @description This function runs pairwise differential CellChat analysis between
#' two biological conditions from one processed, normalized, cell-type-annotated
#' Seurat object. It runs CellChat independently for each condition on shared
#' cell groups, compares inferred communication counts, strengths, pathways,
#' ligand-receptor axes, and sender/receiver roles, then builds a compact
#' LLM-ready differential interpretation context.
#'
#' @author Vijay Nagarajan, PhD, NEI/NIH
#'
#' @param seurat_object_name Character string representing the name of the
#'   Seurat object in the current environment.
#' @param group_by Metadata column containing cell group or cell type labels.
#'   Default is "celltype".
#' @param condition_by Metadata column containing biological condition labels.
#'   Default is "condition".
#' @param condition_a Reference condition for the comparison.
#' @param condition_b Condition compared against `condition_a`.
#' @param species Species for CellChatDB. Supported values are "human" and
#'   "mouse". Default is "human".
#' @param assay Seurat assay containing normalized expression data. Default is
#'   "RNA".
#' @param database_scope CellChatDB scope. Phase 2A supports "default", which
#'   uses the standard CellChatDB subset for protein-mediated communication when
#'   available.
#' @param average_method Average expression method passed to
#'   `CellChat::computeCommunProb`. Default is "triMean".
#' @param trim Optional trim value passed to `CellChat::computeCommunProb`.
#'   Default is NULL.
#' @param raw_use Logical value passed to `CellChat::computeCommunProb`.
#'   Default is TRUE.
#' @param min_cells Minimum cells per group passed to
#'   `CellChat::filterCommunication`. Default is 10.
#' @param pval_threshold P-value threshold for significant CellChat
#'   interactions. Default is 0.05.
#' @param top_n_interactions Number of top differential interactions to include
#'   in the LLM context. Default is 50.
#' @param top_n_pathways Number of top differential pathways to include in the
#'   LLM context. Default is 20.
#' @param run_llm Logical value indicating whether to query an LLM. If FALSE,
#'   no API key is required and `llm_response` is NULL. Default is TRUE.
#' @param output_dir Optional directory for CellChat and summary output files.
#'   Files are only written when this is provided. Default is NULL.
#' @param experimental_context Optional character string describing the
#'   biological experiment for LLM interpretation. Default is NULL.
#' @param llm_server The LLM server to use. Options are "google", "ollama", or
#'   "openai". Default is "google".
#' @param model_G Character string specifying the Gemini model to use for
#'   analysis. Default is "gemini-1.5-flash-latest".
#' @param model_O Character string specifying the Ollama model to use for
#'   analysis. Default is "llama3".
#' @param model_C Character string specifying the OpenAI model to use for
#'   analysis. Default is "gpt-4o-mini".
#' @param api_key_file Character string specifying the path to a file containing
#'   the API key for Google Gemini or OpenAI. Not used when `run_llm = FALSE`
#'   or `llm_server = "ollama"`.
#' @param temperature Numeric value controlling LLM response variability.
#'   Default is 0.
#' @param max_output_tokens Integer specifying the maximum number of tokens the
#'   LLM can generate. Default is 10048.
#' @param model_params A list of parameters passed to `rollama::query` for
#'   Ollama models. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.
#'
#' @return A structured list containing `condition_a`, `condition_b`,
#'   `cellchat_condition_a`, `cellchat_condition_b`, `merged_cellchat`,
#'   `metadata`, per-condition summaries, differential tables, `llm_context`,
#'   `llm_prompt`, and `llm_response`.
#'
#' @usage
#' SCassist_compare_interactions(seurat_object_name,
#'                               group_by = "celltype",
#'                               condition_by = "condition",
#'                               condition_a,
#'                               condition_b,
#'                               species = "human",
#'                               assay = "RNA",
#'                               database_scope = "default",
#'                               average_method = "triMean",
#'                               trim = NULL,
#'                               raw_use = TRUE,
#'                               min_cells = 10,
#'                               pval_threshold = 0.05,
#'                               top_n_interactions = 50,
#'                               top_n_pathways = 20,
#'                               run_llm = TRUE,
#'                               output_dir = NULL,
#'                               experimental_context = NULL,
#'                               llm_server = "google",
#'                               model_G = "gemini-1.5-flash-latest",
#'                               model_O = "llama3",
#'                               model_C = "gpt-4o-mini",
#'                               api_key_file = "api_keys.txt",
#'                               temperature = 0,
#'                               max_output_tokens = 10048,
#'                               model_params = list(seed = 42,
#'                                 temperature = 0, num_gpu = 0))
#'
#' @examples
#' \dontrun{
#' comparison_results <- SCassist_compare_interactions(
#'   seurat_object_name = "seurat_obj",
#'   group_by = "celltype",
#'   condition_by = "condition",
#'   condition_a = "control",
#'   condition_b = "treated",
#'   run_llm = FALSE
#' )
#' }
#'
#' @keywords single-cell, CellChat, differential, ligand-receptor, LLM
#' @export

SCassist_compare_interactions <- function(seurat_object_name,
                                          group_by = "celltype",
                                          condition_by = "condition",
                                          condition_a,
                                          condition_b,
                                          species = "human",
                                          assay = "RNA",
                                          database_scope = "default",
                                          average_method = "triMean",
                                          trim = NULL,
                                          raw_use = TRUE,
                                          min_cells = 10,
                                          pval_threshold = 0.05,
                                          top_n_interactions = 50,
                                          top_n_pathways = 20,
                                          run_llm = TRUE,
                                          output_dir = NULL,
                                          experimental_context = NULL,
                                          llm_server = "google",
                                          model_G = "gemini-1.5-flash-latest",
                                          model_O = "llama3",
                                          model_C = "gpt-4o-mini",
                                          api_key_file = "api_keys.txt",
                                          temperature = 0,
                                          max_output_tokens = 10048,
                                          model_params = list(seed = 42, temperature = 0, num_gpu = 0)
) {

  # input validation
  SCassist_interactions_validate_static_inputs(
    species = species,
    database_scope = database_scope,
    average_method = average_method,
    trim = trim,
    raw_use = raw_use,
    min_cells = min_cells,
    pval_threshold = pval_threshold,
    top_n_interactions = top_n_interactions,
    top_n_pathways = top_n_pathways,
    run_llm = run_llm,
    llm_server = llm_server
  )

  if (missing(condition_a)) {
    stop("Error: condition_a must be provided.", call. = FALSE)
  }
  if (missing(condition_b)) {
    stop("Error: condition_b must be provided.", call. = FALSE)
  }

  condition_args <- SCassist_compare_interactions_validate_condition_args(
    condition_a = condition_a,
    condition_b = condition_b
  )
  condition_a <- condition_args$condition_a
  condition_b <- condition_args$condition_b

  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in environment.", call. = FALSE)
    }
  )

  SCassist_interactions_validate_seurat_inputs(
    seurat_object = seurat_object,
    seurat_object_name = seurat_object_name,
    group_by = group_by,
    assay = assay
  )

  SCassist_compare_interactions_validate_condition_inputs(
    seurat_object = seurat_object,
    condition_by = condition_by,
    condition_a = condition_a,
    condition_b = condition_b
  )

  condition_info <- SCassist_compare_interactions_get_condition_groups(
    seurat_object = seurat_object,
    group_by = group_by,
    condition_by = condition_by,
    condition_a = condition_a,
    condition_b = condition_b
  )

  if (length(condition_info$shared_cell_groups) == 0) {
    stop(
      "Error: no shared cell groups were found between condition_a and condition_b.",
      call. = FALSE
    )
  }

  if (length(condition_info$condition_a_only_cell_groups) > 0 ||
      length(condition_info$condition_b_only_cell_groups) > 0) {
    warning(
      "Phase 2A compares shared cell groups only. Excluded condition-specific groups: ",
      condition_a, " only = ",
      SCassist_compare_interactions_collapse_or_none(condition_info$condition_a_only_cell_groups),
      "; ", condition_b, " only = ",
      SCassist_compare_interactions_collapse_or_none(condition_info$condition_b_only_cell_groups),
      ".",
      call. = FALSE
    )
  }

  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop(
      "Package 'CellChat' is required for SCassist_compare_interactions. ",
      "Please install CellChat v2 from jinworks/CellChat before running this function.",
      call. = FALSE
    )
  }

  condition_objects <- SCassist_compare_interactions_subset_shared_groups(
    seurat_object = seurat_object,
    condition_info = condition_info
  )

  # CellChat execution
  condition_a_results <- SCassist_interactions_run_analysis(
    seurat_object = condition_objects$condition_a_object,
    group_by = group_by,
    species = species,
    assay = assay,
    database_scope = database_scope,
    average_method = average_method,
    trim = trim,
    raw_use = raw_use,
    min_cells = min_cells,
    pval_threshold = pval_threshold,
    top_n_interactions = top_n_interactions,
    top_n_pathways = top_n_pathways,
    experimental_context = experimental_context
  )

  condition_b_results <- SCassist_interactions_run_analysis(
    seurat_object = condition_objects$condition_b_object,
    group_by = group_by,
    species = species,
    assay = assay,
    database_scope = database_scope,
    average_method = average_method,
    trim = trim,
    raw_use = raw_use,
    min_cells = min_cells,
    pval_threshold = pval_threshold,
    top_n_interactions = top_n_interactions,
    top_n_pathways = top_n_pathways,
    experimental_context = experimental_context
  )

  merged_cellchat <- SCassist_compare_interactions_merge_cellchat(
    cellchat_condition_a = condition_a_results$cellchat_object,
    cellchat_condition_b = condition_b_results$cellchat_object,
    condition_a = condition_a,
    condition_b = condition_b
  )

  # summary construction
  metadata <- SCassist_compare_interactions_build_metadata(
    seurat_object = seurat_object,
    condition_info = condition_info,
    condition_objects = condition_objects,
    condition_a_results = condition_a_results,
    condition_b_results = condition_b_results,
    species = species,
    database_scope = database_scope,
    group_by = group_by,
    condition_by = condition_by,
    condition_a = condition_a,
    condition_b = condition_b,
    assay = assay,
    average_method = average_method,
    trim = trim,
    raw_use = raw_use,
    min_cells = min_cells,
    pval_threshold = pval_threshold
  )

  condition_a_summary <- SCassist_compare_interactions_condition_summary(condition_a_results)
  condition_b_summary <- SCassist_compare_interactions_condition_summary(condition_b_results)

  global_comparison <- SCassist_compare_interactions_global_comparison(
    interactions_condition_a = condition_a_results$interactions,
    interactions_condition_b = condition_b_results$interactions,
    condition_a = condition_a,
    condition_b = condition_b
  )

  differential_edges <- SCassist_compare_interactions_compare_edges(
    interactions_condition_a = condition_a_results$interactions,
    interactions_condition_b = condition_b_results$interactions
  )

  differential_pathways <- SCassist_compare_interactions_compare_pathways(
    pathway_summary_condition_a = condition_a_results$pathway_summary,
    pathway_summary_condition_b = condition_b_results$pathway_summary
  )

  differential_interactions <- SCassist_compare_interactions_compare_ligand_receptors(
    interactions_condition_a = condition_a_results$interactions,
    interactions_condition_b = condition_b_results$interactions
  )

  differential_cell_roles <- SCassist_compare_interactions_compare_cell_roles(
    cell_roles_condition_a = condition_a_results$cell_role_summary,
    cell_roles_condition_b = condition_b_results$cell_role_summary
  )

  gained_lost_summary <- SCassist_compare_interactions_gained_lost_summary(
    differential_pathways = differential_pathways,
    differential_interactions = differential_interactions,
    top_n_pathways = top_n_pathways,
    top_n_interactions = top_n_interactions
  )

  # prompt construction
  llm_context <- SCassist_compare_interactions_build_llm_context(
    metadata = metadata,
    global_comparison = global_comparison,
    differential_pathways = differential_pathways,
    differential_interactions = differential_interactions,
    differential_cell_roles = differential_cell_roles,
    gained_lost_summary = gained_lost_summary,
    top_n_interactions = top_n_interactions,
    top_n_pathways = top_n_pathways,
    experimental_context = experimental_context
  )

  llm_prompt <- SCassist_compare_interactions_build_llm_prompt(llm_context)

  comparison_results <- list(
    condition_a = condition_a,
    condition_b = condition_b,
    cellchat_condition_a = condition_a_results$cellchat_object,
    cellchat_condition_b = condition_b_results$cellchat_object,
    merged_cellchat = merged_cellchat,
    metadata = metadata,
    condition_a_summary = condition_a_summary,
    condition_b_summary = condition_b_summary,
    global_comparison = global_comparison,
    differential_edges = differential_edges,
    differential_pathways = differential_pathways,
    differential_interactions = differential_interactions,
    differential_cell_roles = differential_cell_roles,
    gained_lost_summary = gained_lost_summary,
    llm_context = llm_context,
    llm_prompt = llm_prompt,
    llm_response = NULL
  )

  # provider dispatch
  comparison_results <- SCassist_compare_interactions_add_llm_response(
    comparison_results = comparison_results,
    run_llm = run_llm,
    llm_server = llm_server,
    temperature = temperature,
    max_output_tokens = max_output_tokens,
    model_G = model_G,
    model_O = model_O,
    model_C = model_C,
    api_key_file = api_key_file,
    model_params = model_params
  )

  # optional file output
  SCassist_compare_interactions_write_outputs(
    comparison_results = comparison_results,
    output_dir = output_dir,
    run_llm = run_llm
  )

  return(comparison_results)
}


SCassist_compare_interactions_validate_condition_args <- function(condition_a,
                                                                  condition_b) {
  if (!is.character(condition_a) && !is.factor(condition_a)) {
    stop("Error: condition_a must be a single condition label.", call. = FALSE)
  }
  if (!is.character(condition_b) && !is.factor(condition_b)) {
    stop("Error: condition_b must be a single condition label.", call. = FALSE)
  }

  condition_a <- as.character(condition_a)
  condition_b <- as.character(condition_b)

  if (length(condition_a) != 1 || is.na(condition_a) || condition_a == "") {
    stop("Error: condition_a must be a single non-empty condition label.", call. = FALSE)
  }
  if (length(condition_b) != 1 || is.na(condition_b) || condition_b == "") {
    stop("Error: condition_b must be a single non-empty condition label.", call. = FALSE)
  }
  if (identical(condition_a, condition_b)) {
    stop("Error: condition_a and condition_b must be different.", call. = FALSE)
  }

  list(condition_a = condition_a, condition_b = condition_b)
}


SCassist_compare_interactions_validate_condition_inputs <- function(seurat_object,
                                                                    condition_by,
                                                                    condition_a,
                                                                    condition_b) {
  if (!is.character(condition_by) || length(condition_by) != 1 || is.na(condition_by)) {
    stop("Error: condition_by must be a single metadata column name.", call. = FALSE)
  }

  if (!condition_by %in% colnames(seurat_object@meta.data)) {
    stop("Error: condition_by column '", condition_by, "' was not found in seurat_object@meta.data.", call. = FALSE)
  }

  if (any(is.na(seurat_object@meta.data[[condition_by]]))) {
    stop("Error: condition_by column '", condition_by, "' contains NA values.", call. = FALSE)
  }

  conditions <- as.character(seurat_object@meta.data[[condition_by]])
  if (!condition_a %in% conditions) {
    stop("Error: condition_a '", condition_a, "' was not found in the condition_by column.", call. = FALSE)
  }
  if (!condition_b %in% conditions) {
    stop("Error: condition_b '", condition_b, "' was not found in the condition_by column.", call. = FALSE)
  }
}


SCassist_compare_interactions_get_condition_groups <- function(seurat_object,
                                                               group_by,
                                                               condition_by,
                                                               condition_a,
                                                               condition_b) {
  meta <- seurat_object@meta.data
  conditions <- as.character(meta[[condition_by]])
  groups <- as.character(meta[[group_by]])
  cells <- rownames(meta)

  index_a <- which(conditions == condition_a)
  index_b <- which(conditions == condition_b)
  groups_a <- sort(unique(groups[index_a]))
  groups_b <- sort(unique(groups[index_b]))

  shared_groups <- sort(intersect(groups_a, groups_b))
  condition_a_only <- sort(setdiff(groups_a, groups_b))
  condition_b_only <- sort(setdiff(groups_b, groups_a))

  list(
    condition_a_cells_before = cells[index_a],
    condition_b_cells_before = cells[index_b],
    condition_a_cell_groups = groups_a,
    condition_b_cell_groups = groups_b,
    shared_cell_groups = shared_groups,
    condition_a_only_cell_groups = condition_a_only,
    condition_b_only_cell_groups = condition_b_only,
    condition_a_cells_after = cells[conditions == condition_a & groups %in% shared_groups],
    condition_b_cells_after = cells[conditions == condition_b & groups %in% shared_groups],
    condition_a_group_sizes_before = sort(table(groups[index_a]), decreasing = TRUE),
    condition_b_group_sizes_before = sort(table(groups[index_b]), decreasing = TRUE),
    condition_a_group_sizes_after = sort(table(groups[conditions == condition_a & groups %in% shared_groups]), decreasing = TRUE),
    condition_b_group_sizes_after = sort(table(groups[conditions == condition_b & groups %in% shared_groups]), decreasing = TRUE)
  )
}


SCassist_compare_interactions_subset_shared_groups <- function(seurat_object,
                                                               condition_info) {
  condition_a_object <- seurat_object[, condition_info$condition_a_cells_after]
  condition_b_object <- seurat_object[, condition_info$condition_b_cells_after]

  list(
    condition_a_object = condition_a_object,
    condition_b_object = condition_b_object
  )
}


SCassist_compare_interactions_merge_cellchat <- function(cellchat_condition_a,
                                                         cellchat_condition_b,
                                                         condition_a,
                                                         condition_b) {
  tryCatch(
    {
      CellChat::mergeCellChat(
        object.list = stats::setNames(
          list(cellchat_condition_a, cellchat_condition_b),
          c(condition_a, condition_b)
        ),
        add.names = c(condition_a, condition_b)
      )
    },
    error = function(e) {
      warning("CellChat::mergeCellChat failed; returning merged_cellchat = NULL.", call. = FALSE)
      NULL
    }
  )
}


SCassist_compare_interactions_build_metadata <- function(seurat_object,
                                                         condition_info,
                                                         condition_objects,
                                                         condition_a_results,
                                                         condition_b_results,
                                                         species,
                                                         database_scope,
                                                         group_by,
                                                         condition_by,
                                                         condition_a,
                                                         condition_b,
                                                         assay,
                                                         average_method,
                                                         trim,
                                                         raw_use,
                                                         min_cells,
                                                         pval_threshold) {
  database_used <- condition_a_results$metadata$database_used
  database_note <- condition_a_results$metadata$database_note

  list(
    tool = "CellChat",
    analysis_type = "pairwise_condition_comparison",
    cellchat_version = as.character(utils::packageVersion("CellChat")),
    species = species,
    database_used = database_used,
    database_note = database_note,
    database_scope = database_scope,
    group_by = group_by,
    condition_by = condition_by,
    condition_a = condition_a,
    condition_b = condition_b,
    assay = assay,
    average_method = average_method,
    trim = trim,
    raw_use = raw_use,
    min_cells = min_cells,
    pval_threshold = pval_threshold,
    total_cells_original_object = ncol(seurat_object),
    cells_condition_a_before_shared_filtering = length(condition_info$condition_a_cells_before),
    cells_condition_b_before_shared_filtering = length(condition_info$condition_b_cells_before),
    cells_condition_a_after_shared_filtering = ncol(condition_objects$condition_a_object),
    cells_condition_b_after_shared_filtering = ncol(condition_objects$condition_b_object),
    shared_cell_groups = condition_info$shared_cell_groups,
    condition_a_only_cell_groups = condition_info$condition_a_only_cell_groups,
    condition_b_only_cell_groups = condition_info$condition_b_only_cell_groups,
    condition_a_cell_group_sizes_before = as.list(condition_info$condition_a_group_sizes_before),
    condition_b_cell_group_sizes_before = as.list(condition_info$condition_b_group_sizes_before),
    condition_a_cell_group_sizes_after = as.list(condition_info$condition_a_group_sizes_after),
    condition_b_cell_group_sizes_after = as.list(condition_info$condition_b_group_sizes_after),
    generated_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}


SCassist_compare_interactions_condition_summary <- function(interaction_results) {
  list(
    metadata = interaction_results$metadata,
    interactions = interaction_results$interactions,
    pathway_summary = interaction_results$pathway_summary,
    cell_role_summary = interaction_results$cell_role_summary,
    autocrine_paracrine_summary = interaction_results$autocrine_paracrine_summary,
    interaction_type_summary = interaction_results$interaction_type_summary
  )
}


SCassist_compare_interactions_global_comparison <- function(interactions_condition_a,
                                                            interactions_condition_b,
                                                            condition_a,
                                                            condition_b) {
  count_a <- nrow(interactions_condition_a)
  count_b <- nrow(interactions_condition_b)
  strength_a <- SCassist_compare_interactions_sum_score(interactions_condition_a)
  strength_b <- SCassist_compare_interactions_sum_score(interactions_condition_b)
  delta_count <- count_b - count_a
  delta_strength <- strength_b - strength_a

  list(
    total_significant_interactions_condition_a = count_a,
    total_significant_interactions_condition_b = count_b,
    delta_interaction_count = delta_count,
    percent_change_interaction_count = SCassist_compare_interactions_percent_change(count_a, count_b),
    total_communication_strength_condition_a = strength_a,
    total_communication_strength_condition_b = strength_b,
    delta_communication_strength = delta_strength,
    percent_change_communication_strength = SCassist_compare_interactions_percent_change(strength_a, strength_b),
    greater_total_inferred_communication = if (count_b > count_a) condition_b else if (count_a > count_b) condition_a else "equal",
    greater_total_communication_strength = if (strength_b > strength_a) condition_b else if (strength_a > strength_b) condition_a else "equal"
  )
}


SCassist_compare_interactions_compare_edges <- function(interactions_condition_a,
                                                        interactions_condition_b) {
  edges_a <- SCassist_compare_interactions_edge_summary(interactions_condition_a, "condition_a")
  edges_b <- SCassist_compare_interactions_edge_summary(interactions_condition_b, "condition_b")

  edges <- merge(edges_a, edges_b, by = c("source", "target"), all = TRUE)
  if (nrow(edges) == 0) {
    return(SCassist_compare_interactions_empty_edges())
  }

  edges$count_condition_a[is.na(edges$count_condition_a)] <- 0
  edges$count_condition_b[is.na(edges$count_condition_b)] <- 0
  edges$strength_condition_a[is.na(edges$strength_condition_a)] <- 0
  edges$strength_condition_b[is.na(edges$strength_condition_b)] <- 0
  edges$delta_count <- edges$count_condition_b - edges$count_condition_a
  edges$delta_strength <- edges$strength_condition_b - edges$strength_condition_a
  edges$change_direction <- mapply(
    SCassist_compare_interactions_status,
    present_a = edges$count_condition_a > 0,
    present_b = edges$count_condition_b > 0,
    value_a = edges$strength_condition_a,
    value_b = edges$strength_condition_b,
    USE.NAMES = FALSE
  )
  edges$autocrine <- as.character(edges$source) == as.character(edges$target)

  edges <- edges[order(abs(edges$delta_strength), abs(edges$delta_count), decreasing = TRUE), , drop = FALSE]
  rownames(edges) <- NULL
  return(edges)
}


SCassist_compare_interactions_edge_summary <- function(interactions, suffix) {
  count_name <- paste0("count_", suffix)
  strength_name <- paste0("strength_", suffix)

  if (nrow(interactions) == 0) {
    result <- data.frame(source = character(), target = character(), stringsAsFactors = FALSE)
    result[[count_name]] <- numeric()
    result[[strength_name]] <- numeric()
    return(result)
  }

  score <- suppressWarnings(as.numeric(interactions$communication_score))
  score[is.na(score)] <- 0
  edge_data <- data.frame(
    source = as.character(interactions$source),
    target = as.character(interactions$target),
    count = rep(1, nrow(interactions)),
    strength = score,
    stringsAsFactors = FALSE
  )

  result <- aggregate(
    edge_data[, c("count", "strength"), drop = FALSE],
    by = edge_data[, c("source", "target"), drop = FALSE],
    FUN = sum,
    na.rm = TRUE
  )
  names(result)[names(result) == "count"] <- count_name
  names(result)[names(result) == "strength"] <- strength_name
  result
}


SCassist_compare_interactions_empty_edges <- function() {
  data.frame(
    source = character(),
    target = character(),
    count_condition_a = numeric(),
    count_condition_b = numeric(),
    delta_count = numeric(),
    strength_condition_a = numeric(),
    strength_condition_b = numeric(),
    delta_strength = numeric(),
    change_direction = character(),
    autocrine = logical(),
    stringsAsFactors = FALSE
  )
}


SCassist_compare_interactions_compare_pathways <- function(pathway_summary_condition_a,
                                                           pathway_summary_condition_b) {
  pathways_a <- SCassist_compare_interactions_pathway_subset(pathway_summary_condition_a, "condition_a")
  pathways_b <- SCassist_compare_interactions_pathway_subset(pathway_summary_condition_b, "condition_b")

  pathways <- merge(pathways_a, pathways_b, by = "pathway_name", all = TRUE)
  if (nrow(pathways) == 0) {
    return(SCassist_compare_interactions_empty_pathway_comparison())
  }

  numeric_cols <- c(
    "information_flow_condition_a",
    "information_flow_condition_b",
    "n_interactions_condition_a",
    "n_interactions_condition_b"
  )
  for (column in numeric_cols) {
    pathways[[column]][is.na(pathways[[column]])] <- 0
  }

  pathways$delta_information_flow <- pathways$information_flow_condition_b - pathways$information_flow_condition_a
  pathways$percent_change <- mapply(
    SCassist_compare_interactions_percent_change,
    pathways$information_flow_condition_a,
    pathways$information_flow_condition_b,
    USE.NAMES = FALSE
  )
  pathways$status <- mapply(
    SCassist_compare_interactions_status,
    present_a = pathways$n_interactions_condition_a > 0,
    present_b = pathways$n_interactions_condition_b > 0,
    value_a = pathways$information_flow_condition_a,
    value_b = pathways$information_flow_condition_b,
    USE.NAMES = FALSE
  )

  pathways <- pathways[, c(
    "pathway_name",
    "information_flow_condition_a",
    "information_flow_condition_b",
    "delta_information_flow",
    "percent_change",
    "n_interactions_condition_a",
    "n_interactions_condition_b",
    "dominant_senders_condition_a",
    "dominant_senders_condition_b",
    "dominant_receivers_condition_a",
    "dominant_receivers_condition_b",
    "status"
  ), drop = FALSE]

  pathways <- pathways[order(abs(pathways$delta_information_flow), decreasing = TRUE), , drop = FALSE]
  rownames(pathways) <- NULL
  return(pathways)
}


SCassist_compare_interactions_pathway_subset <- function(pathway_summary, suffix) {
  columns <- c(
    "pathway_name",
    "number_of_significant_interactions",
    "total_communication_score",
    "dominant_senders",
    "dominant_receivers"
  )

  if (nrow(pathway_summary) == 0) {
    result <- data.frame(pathway_name = character(), stringsAsFactors = FALSE)
    result[[paste0("information_flow_", suffix)]] <- numeric()
    result[[paste0("n_interactions_", suffix)]] <- numeric()
    result[[paste0("dominant_senders_", suffix)]] <- character()
    result[[paste0("dominant_receivers_", suffix)]] <- character()
    return(result)
  }

  missing_columns <- setdiff(columns, names(pathway_summary))
  for (column in missing_columns) {
    pathway_summary[[column]] <- NA
  }

  result <- data.frame(
    pathway_name = as.character(pathway_summary$pathway_name),
    information_flow = suppressWarnings(as.numeric(pathway_summary$total_communication_score)),
    n_interactions = suppressWarnings(as.numeric(pathway_summary$number_of_significant_interactions)),
    dominant_senders = as.character(pathway_summary$dominant_senders),
    dominant_receivers = as.character(pathway_summary$dominant_receivers),
    stringsAsFactors = FALSE
  )
  result$information_flow[is.na(result$information_flow)] <- 0
  result$n_interactions[is.na(result$n_interactions)] <- 0

  names(result)[names(result) == "information_flow"] <- paste0("information_flow_", suffix)
  names(result)[names(result) == "n_interactions"] <- paste0("n_interactions_", suffix)
  names(result)[names(result) == "dominant_senders"] <- paste0("dominant_senders_", suffix)
  names(result)[names(result) == "dominant_receivers"] <- paste0("dominant_receivers_", suffix)
  result
}


SCassist_compare_interactions_empty_pathway_comparison <- function() {
  data.frame(
    pathway_name = character(),
    information_flow_condition_a = numeric(),
    information_flow_condition_b = numeric(),
    delta_information_flow = numeric(),
    percent_change = numeric(),
    n_interactions_condition_a = numeric(),
    n_interactions_condition_b = numeric(),
    dominant_senders_condition_a = character(),
    dominant_senders_condition_b = character(),
    dominant_receivers_condition_a = character(),
    dominant_receivers_condition_b = character(),
    status = character(),
    stringsAsFactors = FALSE
  )
}


SCassist_compare_interactions_compare_ligand_receptors <- function(interactions_condition_a,
                                                                   interactions_condition_b) {
  interactions_a <- SCassist_compare_interactions_lr_subset(interactions_condition_a, "condition_a")
  interactions_b <- SCassist_compare_interactions_lr_subset(interactions_condition_b, "condition_b")

  interactions <- merge(
    interactions_a,
    interactions_b,
    by = c("source", "target", "interaction_name"),
    all = TRUE
  )

  if (nrow(interactions) == 0) {
    return(SCassist_compare_interactions_empty_lr_comparison())
  }

  present_a <- !is.na(interactions$score_condition_a)
  present_b <- !is.na(interactions$score_condition_b)

  interactions$score_condition_a[is.na(interactions$score_condition_a)] <- 0
  interactions$score_condition_b[is.na(interactions$score_condition_b)] <- 0
  interactions$delta_score <- interactions$score_condition_b - interactions$score_condition_a
  interactions$status <- mapply(
    SCassist_compare_interactions_status,
    present_a = present_a,
    present_b = present_b,
    value_a = interactions$score_condition_a,
    value_b = interactions$score_condition_b,
    USE.NAMES = FALSE
  )

  interactions$ligand <- SCassist_compare_interactions_coalesce_character(
    interactions$ligand_condition_b,
    interactions$ligand_condition_a
  )
  interactions$receptor <- SCassist_compare_interactions_coalesce_character(
    interactions$receptor_condition_b,
    interactions$receptor_condition_a
  )
  interactions$pathway_name <- SCassist_compare_interactions_coalesce_character(
    interactions$pathway_name_condition_b,
    interactions$pathway_name_condition_a
  )
  interactions$annotation <- SCassist_compare_interactions_coalesce_character(
    interactions$annotation_condition_b,
    interactions$annotation_condition_a
  )
  interactions$autocrine <- as.character(interactions$source) == as.character(interactions$target)

  interactions <- interactions[, c(
    "source",
    "target",
    "ligand",
    "receptor",
    "interaction_name",
    "pathway_name",
    "score_condition_a",
    "score_condition_b",
    "delta_score",
    "pval_condition_a",
    "pval_condition_b",
    "status",
    "annotation",
    "autocrine"
  ), drop = FALSE]

  interactions <- interactions[order(abs(interactions$delta_score), decreasing = TRUE), , drop = FALSE]
  rownames(interactions) <- NULL
  return(interactions)
}


SCassist_compare_interactions_lr_subset <- function(interactions, suffix) {
  base_columns <- c(
    "source",
    "target",
    "ligand",
    "receptor",
    "interaction_name",
    "pathway_name",
    "communication_score",
    "pval",
    "annotation"
  )

  if (nrow(interactions) == 0) {
    result <- data.frame(
      source = character(),
      target = character(),
      interaction_name = character(),
      stringsAsFactors = FALSE
    )
    result[[paste0("ligand_", suffix)]] <- character()
    result[[paste0("receptor_", suffix)]] <- character()
    result[[paste0("pathway_name_", suffix)]] <- character()
    result[[paste0("score_", suffix)]] <- numeric()
    result[[paste0("pval_", suffix)]] <- numeric()
    result[[paste0("annotation_", suffix)]] <- character()
    return(result)
  }

  missing_columns <- setdiff(base_columns, names(interactions))
  for (column in missing_columns) {
    interactions[[column]] <- NA
  }

  result <- data.frame(
    source = as.character(interactions$source),
    target = as.character(interactions$target),
    interaction_name = as.character(interactions$interaction_name),
    ligand = as.character(interactions$ligand),
    receptor = as.character(interactions$receptor),
    pathway_name = as.character(interactions$pathway_name),
    score = suppressWarnings(as.numeric(interactions$communication_score)),
    pval = suppressWarnings(as.numeric(interactions$pval)),
    annotation = as.character(interactions$annotation),
    stringsAsFactors = FALSE
  )

  names(result)[names(result) == "ligand"] <- paste0("ligand_", suffix)
  names(result)[names(result) == "receptor"] <- paste0("receptor_", suffix)
  names(result)[names(result) == "pathway_name"] <- paste0("pathway_name_", suffix)
  names(result)[names(result) == "score"] <- paste0("score_", suffix)
  names(result)[names(result) == "pval"] <- paste0("pval_", suffix)
  names(result)[names(result) == "annotation"] <- paste0("annotation_", suffix)
  result
}


SCassist_compare_interactions_empty_lr_comparison <- function() {
  data.frame(
    source = character(),
    target = character(),
    ligand = character(),
    receptor = character(),
    interaction_name = character(),
    pathway_name = character(),
    score_condition_a = numeric(),
    score_condition_b = numeric(),
    delta_score = numeric(),
    pval_condition_a = numeric(),
    pval_condition_b = numeric(),
    status = character(),
    annotation = character(),
    autocrine = logical(),
    stringsAsFactors = FALSE
  )
}


SCassist_compare_interactions_compare_cell_roles <- function(cell_roles_condition_a,
                                                            cell_roles_condition_b) {
  roles_a <- SCassist_compare_interactions_role_subset(cell_roles_condition_a, "condition_a")
  roles_b <- SCassist_compare_interactions_role_subset(cell_roles_condition_b, "condition_b")
  roles <- merge(roles_a, roles_b, by = "cell_group", all = TRUE)

  if (nrow(roles) == 0) {
    return(SCassist_compare_interactions_empty_role_comparison())
  }

  numeric_cols <- c(
    "outgoing_strength_condition_a",
    "outgoing_strength_condition_b",
    "incoming_strength_condition_a",
    "incoming_strength_condition_b"
  )
  for (column in numeric_cols) {
    roles[[column]][is.na(roles[[column]])] <- 0
  }

  roles$delta_outgoing_strength <- roles$outgoing_strength_condition_b - roles$outgoing_strength_condition_a
  roles$delta_incoming_strength <- roles$incoming_strength_condition_b - roles$incoming_strength_condition_a
  roles$role_change_summary <- mapply(
    SCassist_compare_interactions_role_change,
    roles$role_label_condition_a,
    roles$role_label_condition_b,
    roles$delta_outgoing_strength,
    roles$delta_incoming_strength,
    USE.NAMES = FALSE
  )

  roles <- roles[, c(
    "cell_group",
    "outgoing_strength_condition_a",
    "outgoing_strength_condition_b",
    "delta_outgoing_strength",
    "incoming_strength_condition_a",
    "incoming_strength_condition_b",
    "delta_incoming_strength",
    "outdegree_rank_condition_a",
    "outdegree_rank_condition_b",
    "indegree_rank_condition_a",
    "indegree_rank_condition_b",
    "role_label_condition_a",
    "role_label_condition_b",
    "role_change_summary"
  ), drop = FALSE]

  roles <- roles[order(
    abs(roles$delta_outgoing_strength) + abs(roles$delta_incoming_strength),
    decreasing = TRUE
  ), , drop = FALSE]
  rownames(roles) <- NULL
  return(roles)
}


SCassist_compare_interactions_role_subset <- function(cell_roles, suffix) {
  base_columns <- c(
    "cell_group",
    "outgoing_strength",
    "incoming_strength",
    "outdegree_rank",
    "indegree_rank",
    "role_label"
  )

  if (nrow(cell_roles) == 0) {
    result <- data.frame(cell_group = character(), stringsAsFactors = FALSE)
    result[[paste0("outgoing_strength_", suffix)]] <- numeric()
    result[[paste0("incoming_strength_", suffix)]] <- numeric()
    result[[paste0("outdegree_rank_", suffix)]] <- numeric()
    result[[paste0("indegree_rank_", suffix)]] <- numeric()
    result[[paste0("role_label_", suffix)]] <- character()
    return(result)
  }

  missing_columns <- setdiff(base_columns, names(cell_roles))
  for (column in missing_columns) {
    cell_roles[[column]] <- NA
  }

  result <- data.frame(
    cell_group = as.character(cell_roles$cell_group),
    outgoing_strength = suppressWarnings(as.numeric(cell_roles$outgoing_strength)),
    incoming_strength = suppressWarnings(as.numeric(cell_roles$incoming_strength)),
    outdegree_rank = suppressWarnings(as.numeric(cell_roles$outdegree_rank)),
    indegree_rank = suppressWarnings(as.numeric(cell_roles$indegree_rank)),
    role_label = as.character(cell_roles$role_label),
    stringsAsFactors = FALSE
  )

  names(result)[names(result) == "outgoing_strength"] <- paste0("outgoing_strength_", suffix)
  names(result)[names(result) == "incoming_strength"] <- paste0("incoming_strength_", suffix)
  names(result)[names(result) == "outdegree_rank"] <- paste0("outdegree_rank_", suffix)
  names(result)[names(result) == "indegree_rank"] <- paste0("indegree_rank_", suffix)
  names(result)[names(result) == "role_label"] <- paste0("role_label_", suffix)
  result
}


SCassist_compare_interactions_empty_role_comparison <- function() {
  data.frame(
    cell_group = character(),
    outgoing_strength_condition_a = numeric(),
    outgoing_strength_condition_b = numeric(),
    delta_outgoing_strength = numeric(),
    incoming_strength_condition_a = numeric(),
    incoming_strength_condition_b = numeric(),
    delta_incoming_strength = numeric(),
    outdegree_rank_condition_a = numeric(),
    outdegree_rank_condition_b = numeric(),
    indegree_rank_condition_a = numeric(),
    indegree_rank_condition_b = numeric(),
    role_label_condition_a = character(),
    role_label_condition_b = character(),
    role_change_summary = character(),
    stringsAsFactors = FALSE
  )
}


SCassist_compare_interactions_gained_lost_summary <- function(differential_pathways,
                                                              differential_interactions,
                                                              top_n_pathways,
                                                              top_n_interactions) {
  records <- list(
    SCassist_compare_interactions_summary_records(
      differential_pathways,
      record_type = "pathway",
      summary_category = "top_gained_pathways_in_condition_b",
      status = "gained_in_condition_b",
      order_column = "information_flow_condition_b",
      decreasing = TRUE,
      n = top_n_pathways
    ),
    SCassist_compare_interactions_summary_records(
      differential_pathways,
      record_type = "pathway",
      summary_category = "top_lost_pathways_in_condition_b",
      status = "lost_in_condition_b",
      order_column = "information_flow_condition_a",
      decreasing = TRUE,
      n = top_n_pathways
    ),
    SCassist_compare_interactions_summary_records(
      differential_pathways,
      record_type = "pathway",
      summary_category = "top_increased_pathways_in_condition_b",
      status = "increased_in_condition_b",
      order_column = "delta_information_flow",
      decreasing = TRUE,
      n = top_n_pathways
    ),
    SCassist_compare_interactions_summary_records(
      differential_pathways,
      record_type = "pathway",
      summary_category = "top_decreased_pathways_in_condition_b",
      status = "decreased_in_condition_b",
      order_column = "delta_information_flow",
      decreasing = FALSE,
      n = top_n_pathways
    ),
    SCassist_compare_interactions_summary_records(
      differential_interactions,
      record_type = "ligand_receptor",
      summary_category = "top_gained_ligand_receptor_interactions_in_condition_b",
      status = "gained_in_condition_b",
      order_column = "score_condition_b",
      decreasing = TRUE,
      n = top_n_interactions
    ),
    SCassist_compare_interactions_summary_records(
      differential_interactions,
      record_type = "ligand_receptor",
      summary_category = "top_lost_ligand_receptor_interactions_in_condition_b",
      status = "lost_in_condition_b",
      order_column = "score_condition_a",
      decreasing = TRUE,
      n = top_n_interactions
    ),
    SCassist_compare_interactions_summary_records(
      differential_interactions,
      record_type = "ligand_receptor",
      summary_category = "top_increased_ligand_receptor_interactions_in_condition_b",
      status = "increased_in_condition_b",
      order_column = "delta_score",
      decreasing = TRUE,
      n = top_n_interactions
    ),
    SCassist_compare_interactions_summary_records(
      differential_interactions,
      record_type = "ligand_receptor",
      summary_category = "top_decreased_ligand_receptor_interactions_in_condition_b",
      status = "decreased_in_condition_b",
      order_column = "delta_score",
      decreasing = FALSE,
      n = top_n_interactions
    )
  )

  records <- records[vapply(records, nrow, integer(1)) > 0]
  if (length(records) == 0) {
    return(SCassist_compare_interactions_empty_gained_lost_summary())
  }

  result <- do.call(rbind, records)
  rownames(result) <- NULL
  return(result)
}


SCassist_compare_interactions_summary_records <- function(data,
                                                          record_type,
                                                          summary_category,
                                                          status,
                                                          order_column,
                                                          decreasing,
                                                          n) {
  if (nrow(data) == 0 || !status %in% data$status || !order_column %in% names(data)) {
    return(SCassist_compare_interactions_empty_gained_lost_summary())
  }

  rows <- data[data$status == status, , drop = FALSE]
  rows <- rows[order(rows[[order_column]], decreasing = decreasing, na.last = TRUE), , drop = FALSE]
  rows <- head(rows, n)

  if (identical(record_type, "pathway")) {
    return(data.frame(
      summary_category = summary_category,
      record_type = record_type,
      pathway_name = rows$pathway_name,
      interaction_name = NA_character_,
      source = NA_character_,
      target = NA_character_,
      ligand = NA_character_,
      receptor = NA_character_,
      score_condition_a = NA_real_,
      score_condition_b = NA_real_,
      delta_score = NA_real_,
      information_flow_condition_a = rows$information_flow_condition_a,
      information_flow_condition_b = rows$information_flow_condition_b,
      delta_information_flow = rows$delta_information_flow,
      status = rows$status,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    summary_category = summary_category,
    record_type = record_type,
    pathway_name = rows$pathway_name,
    interaction_name = rows$interaction_name,
    source = rows$source,
    target = rows$target,
    ligand = rows$ligand,
    receptor = rows$receptor,
    score_condition_a = rows$score_condition_a,
    score_condition_b = rows$score_condition_b,
    delta_score = rows$delta_score,
    information_flow_condition_a = NA_real_,
    information_flow_condition_b = NA_real_,
    delta_information_flow = NA_real_,
    status = rows$status,
    stringsAsFactors = FALSE
  )
}


SCassist_compare_interactions_empty_gained_lost_summary <- function() {
  data.frame(
    summary_category = character(),
    record_type = character(),
    pathway_name = character(),
    interaction_name = character(),
    source = character(),
    target = character(),
    ligand = character(),
    receptor = character(),
    score_condition_a = numeric(),
    score_condition_b = numeric(),
    delta_score = numeric(),
    information_flow_condition_a = numeric(),
    information_flow_condition_b = numeric(),
    delta_information_flow = numeric(),
    status = character(),
    stringsAsFactors = FALSE
  )
}


SCassist_compare_interactions_build_llm_context <- function(metadata,
                                                           global_comparison,
                                                           differential_pathways,
                                                           differential_interactions,
                                                           differential_cell_roles,
                                                           gained_lost_summary,
                                                           top_n_interactions,
                                                           top_n_pathways,
                                                           experimental_context = NULL) {
  list(
    metadata = metadata,
    global_comparison = global_comparison,
    shared_and_missing_cell_groups = list(
      shared_cell_groups = metadata$shared_cell_groups,
      condition_a_only_cell_groups = metadata$condition_a_only_cell_groups,
      condition_b_only_cell_groups = metadata$condition_b_only_cell_groups
    ),
    top_increased_decreased_pathways = head(
      differential_pathways[differential_pathways$status %in% c("increased_in_condition_b", "decreased_in_condition_b"), , drop = FALSE],
      top_n_pathways
    ),
    top_gained_lost_pathways = head(
      differential_pathways[differential_pathways$status %in% c("gained_in_condition_b", "lost_in_condition_b"), , drop = FALSE],
      top_n_pathways
    ),
    top_increased_decreased_ligand_receptor_interactions = head(
      differential_interactions[differential_interactions$status %in% c("increased_in_condition_b", "decreased_in_condition_b"), , drop = FALSE],
      top_n_interactions
    ),
    top_gained_lost_ligand_receptor_interactions = head(
      differential_interactions[differential_interactions$status %in% c("gained_in_condition_b", "lost_in_condition_b"), , drop = FALSE],
      top_n_interactions
    ),
    major_sender_receiver_shifts = SCassist_compare_interactions_major_role_shifts(
      differential_cell_roles = differential_cell_roles,
      n = top_n_pathways
    ),
    gained_lost_summary = gained_lost_summary,
    experimental_context = experimental_context,
    interpretation_caveats = c(
      "CellChat infers potential communication from mRNA expression and curated prior knowledge.",
      "CellChat communication probability is an interaction strength score, not a literal probability.",
      "Differences between conditions do not prove causal signaling changes.",
      "Protein-level activity and spatial proximity are not directly measured from standard scRNA-seq.",
      "Differences may be influenced by cell group definitions, sequencing depth, cell numbers, normalization, and average expression method.",
      "Phase 2A compares shared cell groups only; condition-specific cell groups are reported but not included in the main pairwise communication comparison.",
      "Differential interpretation is based on CellChat communication score / information-flow differences, not DEG-supported ligand/receptor logFC analysis."
    )
  )
}


SCassist_compare_interactions_major_role_shifts <- function(differential_cell_roles, n = 10) {
  if (nrow(differential_cell_roles) == 0) {
    return(differential_cell_roles)
  }

  score <- abs(differential_cell_roles$delta_outgoing_strength) +
    abs(differential_cell_roles$delta_incoming_strength)
  differential_cell_roles <- differential_cell_roles[order(score, decreasing = TRUE), , drop = FALSE]
  head(differential_cell_roles, n)
}


SCassist_compare_interactions_build_llm_prompt <- function(llm_context) {
  context_text <- SCassist_interactions_context_to_text(llm_context)

  # ---- Differential InteractionAgent LLM Prompt ----
  prompt_text <- paste0(
    "You are SCassist InteractionAgent, an AI assistant that interprets CellChat v2 pairwise differential cell-cell communication results for scientists.\n\n",
    "Using only the supplied differential CellChat context, write a concise biological interpretation comparing condition_b against condition_a. Focus on the strongest and most biologically relevant changes rather than listing every result.\n\n",
    "Use cautious scientific language. Prefer phrases such as “CellChat inferred,” “the data suggest,” and “this is consistent with.” Avoid definitive or causal language such as “this proves,” “these cells definitely communicate,” or “this pathway causes.”\n\n",
    "Your response should include:\n\n",
    "1. Overall comparison\n",
    "   - Summarize whether condition_b shows increased, decreased, or redistributed inferred communication compared with condition_a.\n",
    "   - Mention changes in total interaction count and total communication strength.\n\n",
    "2. Major pathway changes\n",
    "   - Highlight pathways that are gained, lost, increased, or decreased in condition_b.\n",
    "   - Use pathway information flow as the main pathway-level metric.\n\n",
    "3. Key ligand-receptor axes\n",
    "   - Identify the strongest differential ligand-receptor pairs and the sender-receiver cell groups they connect.\n",
    "   - Separate gained/lost interactions from increased/decreased shared interactions.\n\n",
    "4. Sender and receiver shifts\n",
    "   - Describe which cell groups show increased or decreased outgoing signaling.\n",
    "   - Describe which cell groups show increased or decreased incoming signaling.\n\n",
    "5. Biological interpretation\n",
    "   - Explain what biological processes the strongest differential pathways may suggest, using the experimental context when available.\n",
    "   - Do not infer mechanisms that are not supported by the supplied context.\n\n",
    "6. Confidence and caveats\n",
    "   - CellChat infers potential communication from mRNA expression and curated ligand-receptor knowledge.\n",
    "   - Communication probability is an interaction strength score, not a literal probability.\n",
    "   - Differences between conditions do not prove causality, protein-level activity, or spatial proximity.\n",
    "   - Gained or lost interactions reflect differences in CellChat-detected/significant communication under the chosen thresholds, not absolute biological presence or absence.\n",
    "   - Phase 2A compares shared cell groups only; condition-specific cell groups are reported but excluded from the main pairwise comparison.\n",
    "   - Differential interpretation is based on CellChat communication scores and pathway information flow, not DEG-supported ligand/receptor logFC analysis.\n\n",
    "7. Suggested follow-up\n",
    "   - Recommend practical validation or follow-up analyses, such as checking ligand/receptor expression, validating key pathways experimentally, comparing additional conditions, or running DEG-supported ligand-receptor analysis.\n\n",
    "Do not invent unsupported pathways, cell types, ligand-receptor pairs, or mechanisms."
  )

  paste0(
    prompt_text,
    "\n\n",
    "Structured differential CellChat context:\n",
    context_text
  )
}


SCassist_compare_interactions_add_llm_response <- function(comparison_results,
                                                          run_llm,
                                                          llm_server,
                                                          temperature,
                                                          max_output_tokens,
                                                          model_G,
                                                          model_O,
                                                          model_C,
                                                          api_key_file,
                                                          model_params) {
  if (!run_llm) {
    comparison_results["llm_response"] <- list(NULL)
    return(comparison_results)
  }

  if (llm_server == "google") {
    comparison_results$llm_response <- SCassist_compare_interactions_G(
      llm_prompt = comparison_results$llm_prompt,
      temperature = temperature,
      max_output_tokens = max_output_tokens,
      model = model_G,
      api_key_file = api_key_file
    )
  } else if (llm_server == "ollama") {
    comparison_results$llm_response <- SCassist_compare_interactions_L(
      llm_prompt = comparison_results$llm_prompt,
      model_params = model_params,
      model = model_O
    )
  } else if (llm_server == "openai") {
    comparison_results$llm_response <- SCassist_compare_interactions_C(
      llm_prompt = comparison_results$llm_prompt,
      temperature = temperature,
      max_output_tokens = max_output_tokens,
      model = model_C,
      api_key_file = api_key_file
    )
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama' or 'openai'.", call. = FALSE)
  }

  return(comparison_results)
}


SCassist_compare_interactions_G <- function(llm_prompt,
                                            temperature = 0,
                                            max_output_tokens = 10048,
                                            model = "gemini-1.5-flash-latest",
                                            api_key_file = "api_keys.txt") {
  SCassist_analyze_interactions_G(
    llm_prompt = llm_prompt,
    temperature = temperature,
    max_output_tokens = max_output_tokens,
    model = model,
    api_key_file = api_key_file
  )
}


SCassist_compare_interactions_L <- function(llm_prompt,
                                            model_params = list(seed = 42, temperature = 0, num_gpu = 0),
                                            model = "llama3") {
  SCassist_analyze_interactions_L(
    llm_prompt = llm_prompt,
    model_params = model_params,
    model = model
  )
}


SCassist_compare_interactions_C <- function(llm_prompt,
                                            temperature = 0,
                                            max_output_tokens = 10048,
                                            model = "gpt-4o-mini",
                                            api_key_file = "api_keys.txt") {
  SCassist_analyze_interactions_C(
    llm_prompt = llm_prompt,
    temperature = temperature,
    max_output_tokens = max_output_tokens,
    model = model,
    api_key_file = api_key_file
  )
}


SCassist_compare_interactions_write_outputs <- function(comparison_results,
                                                       output_dir,
                                                       run_llm) {
  if (is.null(output_dir)) {
    return(invisible(NULL))
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  saveRDS(comparison_results$cellchat_condition_a, file = file.path(output_dir, "cellchat_condition_a.rds"))
  saveRDS(comparison_results$cellchat_condition_b, file = file.path(output_dir, "cellchat_condition_b.rds"))
  saveRDS(comparison_results$merged_cellchat, file = file.path(output_dir, "merged_cellchat.rds"))
  utils::write.csv(comparison_results$condition_a_summary$interactions, file = file.path(output_dir, "condition_a_interactions.csv"), row.names = FALSE)
  utils::write.csv(comparison_results$condition_b_summary$interactions, file = file.path(output_dir, "condition_b_interactions.csv"), row.names = FALSE)
  utils::write.csv(comparison_results$condition_a_summary$pathway_summary, file = file.path(output_dir, "condition_a_pathway_summary.csv"), row.names = FALSE)
  utils::write.csv(comparison_results$condition_b_summary$pathway_summary, file = file.path(output_dir, "condition_b_pathway_summary.csv"), row.names = FALSE)
  utils::write.csv(comparison_results$condition_a_summary$cell_role_summary, file = file.path(output_dir, "condition_a_cell_role_summary.csv"), row.names = FALSE)
  utils::write.csv(comparison_results$condition_b_summary$cell_role_summary, file = file.path(output_dir, "condition_b_cell_role_summary.csv"), row.names = FALSE)
  utils::write.csv(as.data.frame(comparison_results$global_comparison, stringsAsFactors = FALSE), file = file.path(output_dir, "differential_global_comparison.csv"), row.names = FALSE)
  utils::write.csv(comparison_results$differential_edges, file = file.path(output_dir, "differential_edges.csv"), row.names = FALSE)
  utils::write.csv(comparison_results$differential_pathways, file = file.path(output_dir, "differential_pathways.csv"), row.names = FALSE)
  utils::write.csv(comparison_results$differential_interactions, file = file.path(output_dir, "differential_interactions.csv"), row.names = FALSE)
  utils::write.csv(comparison_results$differential_cell_roles, file = file.path(output_dir, "differential_cell_roles.csv"), row.names = FALSE)
  utils::write.csv(comparison_results$gained_lost_summary, file = file.path(output_dir, "gained_lost_summary.csv"), row.names = FALSE)

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(
      comparison_results$llm_context,
      path = file.path(output_dir, "differential_llm_context.json"),
      pretty = TRUE,
      auto_unbox = TRUE,
      na = "null"
    )
  } else {
    warning("Package 'jsonlite' is not available; skipping differential_llm_context.json output.")
  }

  if (run_llm && !is.null(comparison_results$llm_response)) {
    writeLines(comparison_results$llm_response, con = file.path(output_dir, "differential_interpretation.txt"))
  }

  return(invisible(NULL))
}


SCassist_compare_interactions_sum_score <- function(interactions) {
  if (nrow(interactions) == 0 || !"communication_score" %in% names(interactions)) {
    return(0)
  }

  score <- suppressWarnings(as.numeric(interactions$communication_score))
  sum(score, na.rm = TRUE)
}


SCassist_compare_interactions_percent_change <- function(value_a, value_b) {
  value_a <- suppressWarnings(as.numeric(value_a))
  value_b <- suppressWarnings(as.numeric(value_b))

  if (length(value_a) != 1 || length(value_b) != 1 ||
      is.na(value_a) || is.na(value_b) || value_a == 0) {
    return(NA_real_)
  }

  100 * (value_b - value_a) / abs(value_a)
}


SCassist_compare_interactions_status <- function(present_a,
                                                 present_b,
                                                 value_a,
                                                 value_b,
                                                 tolerance = sqrt(.Machine$double.eps)) {
  if (!isTRUE(present_a) && isTRUE(present_b)) {
    return("gained_in_condition_b")
  }
  if (isTRUE(present_a) && !isTRUE(present_b)) {
    return("lost_in_condition_b")
  }

  delta <- suppressWarnings(as.numeric(value_b) - as.numeric(value_a))
  if (is.na(delta) || abs(delta) <= tolerance) {
    return("maintained")
  }
  if (delta > 0) {
    return("increased_in_condition_b")
  }
  return("decreased_in_condition_b")
}


SCassist_compare_interactions_coalesce_character <- function(primary, fallback) {
  primary <- as.character(primary)
  fallback <- as.character(fallback)
  missing_primary <- is.na(primary) | primary == ""
  primary[missing_primary] <- fallback[missing_primary]
  primary
}


SCassist_compare_interactions_role_change <- function(role_a,
                                                      role_b,
                                                      delta_outgoing,
                                                      delta_incoming) {
  outgoing_label <- if (is.na(delta_outgoing) || abs(delta_outgoing) <= sqrt(.Machine$double.eps)) {
    "outgoing_maintained"
  } else if (delta_outgoing > 0) {
    "outgoing_increased_in_condition_b"
  } else {
    "outgoing_decreased_in_condition_b"
  }

  incoming_label <- if (is.na(delta_incoming) || abs(delta_incoming) <= sqrt(.Machine$double.eps)) {
    "incoming_maintained"
  } else if (delta_incoming > 0) {
    "incoming_increased_in_condition_b"
  } else {
    "incoming_decreased_in_condition_b"
  }

  role_label <- if (!is.na(role_a) && !is.na(role_b) && identical(role_a, role_b)) {
    paste0("role_maintained_", role_b)
  } else {
    paste0("role_changed_from_", role_a, "_to_", role_b)
  }

  paste(role_label, outgoing_label, incoming_label, sep = "; ")
}


SCassist_compare_interactions_collapse_or_none <- function(values) {
  if (length(values) == 0) {
    return("none")
  }
  paste(values, collapse = ", ")
}
