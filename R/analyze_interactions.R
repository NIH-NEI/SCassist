#' @title Analyze Cell-Cell Communication with CellChat and an Optional LLM
#'
#' @description This function runs a single-condition CellChat analysis on a
#' processed, normalized, cell-type-annotated Seurat object. It extracts
#' ligand-receptor interactions, pathway-level communication summaries,
#' sender/receiver cell role summaries, autocrine versus paracrine signaling,
#' and CellChatDB interaction annotations where available. It then builds a
#' compact LLM-ready context and can optionally query Gemini, Ollama, or OpenAI
#' using the same provider pattern used by other SCassist functions.
#'
#' @author Vijay Nagarajan, PhD, NEI/NIH
#'
#' @param seurat_object_name Character string representing the name of the
#'   Seurat object in the current environment.
#' @param group_by Metadata column containing cell group or cell type labels.
#'   Default is "celltype".
#' @param species Species for CellChatDB. Supported values are "human" and
#'   "mouse". Default is "human".
#' @param assay Seurat assay containing normalized expression data. Default is
#'   "RNA".
#' @param database_scope CellChatDB scope. Phase 1 supports "default", which
#'   uses protein-mediated CellChatDB interactions when annotations are
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
#' @param top_n_interactions Number of top interactions to include in the LLM
#'   context. Default is 50.
#' @param top_n_pathways Number of top pathways to include in the LLM context.
#'   Default is 20.
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
#' @return A structured list containing `cellchat_object`, `metadata`,
#'   `interactions`, `pathway_summary`, `cell_role_summary`,
#'   `autocrine_paracrine_summary`, `interaction_type_summary`, `llm_context`,
#'   `llm_prompt`, and `llm_response`.
#'
#' @usage
#' SCassist_analyze_interactions(seurat_object_name,
#'                               group_by = "celltype",
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
#' interaction_results <- SCassist_analyze_interactions(
#'   seurat_object_name = "seurat_obj",
#'   group_by = "celltype",
#'   species = "human",
#'   experimental_context = "PBMC single-cell RNA-seq after stimulation",
#'   api_key_file = "api_keys.txt"
#' )
#'
#' interaction_results_no_llm <- SCassist_analyze_interactions(
#'   seurat_object_name = "seurat_obj",
#'   group_by = "celltype",
#'   run_llm = FALSE
#' )
#' }
#'
#' @keywords single-cell, CellChat, ligand-receptor, cell-cell communication, LLM
#' @export

SCassist_analyze_interactions <- function(seurat_object_name,
                                          group_by = "celltype",
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

  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop(
      "Package 'CellChat' is required for SCassist_analyze_interactions. ",
      "Please install CellChat v2 from jinworks/CellChat before running this function.",
      call. = FALSE
    )
  }

  # CellChat execution and structured result extraction
  interaction_results <- SCassist_interactions_run_analysis(
    seurat_object = seurat_object,
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

  # provider dispatch
  interaction_results <- SCassist_interactions_add_llm_response(
    interaction_results = interaction_results,
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
  SCassist_interactions_write_outputs(
    interaction_results = interaction_results,
    output_dir = output_dir,
    run_llm = run_llm
  )

  return(interaction_results)
}


SCassist_interactions_validate_static_inputs <- function(species,
                                                         database_scope,
                                                         average_method,
                                                         trim,
                                                         raw_use,
                                                         min_cells,
                                                         pval_threshold,
                                                         top_n_interactions,
                                                         top_n_pathways,
                                                         run_llm,
                                                         llm_server) {
  valid_species <- c("human", "mouse")
  if (!is.character(species) || length(species) != 1 || is.na(species) ||
      !species %in% valid_species) {
    stop("Error: species must be one of 'human' or 'mouse'.", call. = FALSE)
  }

  if (!is.character(database_scope) || length(database_scope) != 1 ||
      is.na(database_scope) || !database_scope %in% c("default")) {
    stop("Error: unsupported database_scope. Phase 1 supports only 'default'.", call. = FALSE)
  }

  valid_average_methods <- c("triMean", "truncatedMean", "thresholdedMean", "median")
  if (!is.character(average_method) || length(average_method) != 1 ||
      is.na(average_method) || !average_method %in% valid_average_methods) {
    stop(
      "Error: average_method must be one of 'triMean', 'truncatedMean', ",
      "'thresholdedMean', or 'median'.",
      call. = FALSE
    )
  }

  if (!is.null(trim) && (!is.numeric(trim) || length(trim) != 1 || is.na(trim) || trim < 0 || trim > 1)) {
    stop("Error: trim must be NULL or a numeric value between 0 and 1.", call. = FALSE)
  }

  if (!is.logical(raw_use) || length(raw_use) != 1 || is.na(raw_use)) {
    stop("Error: raw_use must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(min_cells) || length(min_cells) != 1 || is.na(min_cells) || min_cells <= 0) {
    stop("Error: min_cells must be a positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(pval_threshold) || length(pval_threshold) != 1 ||
      is.na(pval_threshold) || pval_threshold < 0 || pval_threshold > 1) {
    stop("Error: pval_threshold must be a numeric value between 0 and 1.", call. = FALSE)
  }

  if (!is.numeric(top_n_interactions) || length(top_n_interactions) != 1 ||
      is.na(top_n_interactions) || top_n_interactions <= 0) {
    stop("Error: top_n_interactions must be a positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(top_n_pathways) || length(top_n_pathways) != 1 ||
      is.na(top_n_pathways) || top_n_pathways <= 0) {
    stop("Error: top_n_pathways must be a positive numeric value.", call. = FALSE)
  }

  if (!is.logical(run_llm) || length(run_llm) != 1 || is.na(run_llm)) {
    stop("Error: run_llm must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.character(llm_server) || length(llm_server) != 1 ||
      is.na(llm_server) || !llm_server %in% c("google", "ollama", "openai")) {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama' or 'openai'.", call. = FALSE)
  }
}


SCassist_interactions_validate_seurat_inputs <- function(seurat_object,
                                                         seurat_object_name,
                                                         group_by,
                                                         assay) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for this function. Please install it.", call. = FALSE)
  }

  if (!inherits(seurat_object, "Seurat")) {
    stop("Error: object '", seurat_object_name, "' is not a Seurat object.", call. = FALSE)
  }

  if (!is.character(group_by) || length(group_by) != 1 || is.na(group_by)) {
    stop("Error: group_by must be a single metadata column name.", call. = FALSE)
  }

  if (!group_by %in% colnames(seurat_object@meta.data)) {
    stop("Error: group_by column '", group_by, "' was not found in seurat_object@meta.data.", call. = FALSE)
  }

  if (any(is.na(seurat_object@meta.data[[group_by]]))) {
    stop("Error: group_by column '", group_by, "' contains NA values.", call. = FALSE)
  }

  if (!is.character(assay) || length(assay) != 1 || is.na(assay)) {
    stop("Error: assay must be a single assay name.", call. = FALSE)
  }

  if (!assay %in% names(seurat_object@assays)) {
    stop("Error: assay '", assay, "' was not found in the Seurat object.", call. = FALSE)
  }
}


SCassist_interactions_run_analysis <- function(seurat_object,
                                               group_by,
                                               species,
                                               assay,
                                               database_scope,
                                               average_method,
                                               trim,
                                               raw_use,
                                               min_cells,
                                               pval_threshold,
                                               top_n_interactions,
                                               top_n_pathways,
                                               experimental_context) {
  # CellChat execution
  database_info <- SCassist_interactions_get_cellchat_db(
    species = species,
    database_scope = database_scope
  )

  expression_data <- SCassist_interactions_get_expression_data(
    seurat_object = seurat_object,
    assay = assay
  )

  meta <- seurat_object@meta.data
  meta$SCassist_CellChat_group <- as.character(meta[[group_by]])

  cellchat <- CellChat::createCellChat(
    object = expression_data,
    meta = meta,
    group.by = "SCassist_CellChat_group"
  )

  cellchat@DB <- database_info$database
  cellchat <- CellChat::subsetData(cellchat)
  # Avoid requiring CellChat's optional presto dependency for the fast Wilcoxon path.
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)

  compute_args <- list(
    object = cellchat,
    type = average_method,
    raw.use = raw_use
  )

  if (!is.null(trim)) {
    compute_args$trim <- trim
  }

  cellchat <- do.call(CellChat::computeCommunProb, compute_args)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  cellchat <- SCassist_interactions_compute_centrality(cellchat)

  # result extraction
  interactions <- SCassist_interactions_extract_communications(
    cellchat = cellchat,
    cellchat_db = database_info$database,
    pval_threshold = pval_threshold
  )

  cell_group_sizes <- sort(table(seurat_object@meta.data[[group_by]]), decreasing = TRUE)
  metadata <- SCassist_interactions_build_metadata(
    species = species,
    database_info = database_info,
    group_by = group_by,
    assay = assay,
    average_method = average_method,
    trim = trim,
    raw_use = raw_use,
    min_cells = min_cells,
    pval_threshold = pval_threshold,
    number_of_cells = ncol(seurat_object),
    cell_group_sizes = cell_group_sizes
  )

  # summary construction
  pathway_summary <- SCassist_interactions_summarize_pathways(
    interactions = interactions,
    cell_groups = names(cell_group_sizes)
  )

  cell_role_summary <- SCassist_interactions_summarize_roles(
    interactions = interactions,
    cell_groups = names(cell_group_sizes),
    cellchat = cellchat
  )

  autocrine_paracrine_summary <- SCassist_interactions_summarize_autocrine_paracrine(
    interactions = interactions
  )

  interaction_type_summary <- SCassist_interactions_summarize_interaction_types(
    interactions = interactions
  )

  # prompt construction
  llm_context <- SCassist_interactions_build_llm_context(
    metadata = metadata,
    interactions = interactions,
    pathway_summary = pathway_summary,
    cell_role_summary = cell_role_summary,
    autocrine_paracrine_summary = autocrine_paracrine_summary,
    interaction_type_summary = interaction_type_summary,
    top_n_interactions = top_n_interactions,
    top_n_pathways = top_n_pathways,
    experimental_context = experimental_context
  )

  llm_prompt <- SCassist_interactions_build_llm_prompt(llm_context)

  list(
    cellchat_object = cellchat,
    metadata = metadata,
    interactions = interactions,
    pathway_summary = pathway_summary,
    cell_role_summary = cell_role_summary,
    autocrine_paracrine_summary = autocrine_paracrine_summary,
    interaction_type_summary = interaction_type_summary,
    llm_context = llm_context,
    llm_prompt = llm_prompt,
    llm_response = NULL
  )
}


SCassist_interactions_get_expression_data <- function(seurat_object, assay) {
  expression_data <- tryCatch(
    {
      Seurat::GetAssayData(seurat_object, assay = assay, slot = "data")
    },
    error = function(e) {
      tryCatch(
        {
          Seurat::GetAssayData(seurat_object, assay = assay, layer = "data")
        },
        error = function(e2) {
          NULL
        }
      )
    }
  )

  if (is.null(expression_data) || length(expression_data) == 0) {
    expression_data <- tryCatch(
      {
        Seurat::GetAssayData(seurat_object, assay = assay, slot = "counts")
      },
      error = function(e) {
        tryCatch(
          {
            Seurat::GetAssayData(seurat_object, assay = assay, layer = "counts")
          },
          error = function(e2) {
            NULL
          }
        )
      }
    )
  }

  if (is.null(expression_data) || length(expression_data) == 0) {
    stop("Error: unable to extract expression data from assay '", assay, "'.", call. = FALSE)
  }

  return(expression_data)
}


SCassist_interactions_get_cellchat_db <- function(species, database_scope) {
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop(
      "Package 'CellChat' is required to load CellChatDB. ",
      "Please install CellChat v2 from jinworks/CellChat before running this function.",
      call. = FALSE
    )
  }

  db_name <- if (species == "human") {
    "CellChatDB.human"
  } else {
    "CellChatDB.mouse"
  }

  db_env <- new.env(parent = emptyenv())
  loaded <- tryCatch(
    {
      suppressWarnings(utils::data(list = db_name, package = "CellChat", envir = db_env))
      TRUE
    },
    error = function(e) {
      FALSE
    }
  )

  if (!loaded || !exists(db_name, envir = db_env, inherits = FALSE)) {
    stop("Could not load ", db_name, " from CellChat.", call. = FALSE)
  }

  cellchat_db <- get(db_name, envir = db_env, inherits = FALSE)
  database_note <- "default CellChatDB"

  if (database_scope == "default") {
    db_subset <- SCassist_interactions_subset_default_database(cellchat_db)
    cellchat_db <- db_subset$database
    database_note <- db_subset$note
  }

  list(
    database = cellchat_db,
    database_name = db_name,
    database_note = database_note
  )
}


SCassist_interactions_subset_default_database <- function(cellchat_db) {
  subset_db <- tryCatch(
    {
      CellChat::subsetDB(cellchat_db)
    },
    error = function(e) {
      NULL
    }
  )

  if (!is.null(subset_db) && !is.null(subset_db$interaction) &&
      nrow(subset_db$interaction) > 0) {
    return(list(
      database = subset_db,
      note = "default CellChatDB subset generated with CellChat::subsetDB"
    ))
  }

  note <- "default protein-mediated CellChatDB interactions"

  if (is.null(cellchat_db$interaction) || is.null(cellchat_db$interaction$annotation)) {
    return(list(
      database = cellchat_db,
      note = "default CellChatDB; interaction annotation column was unavailable for protein-mediated subsetting"
    ))
  }

  protein_categories <- c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
  keep <- cellchat_db$interaction$annotation %in% protein_categories

  if (!any(keep)) {
    return(list(
      database = cellchat_db,
      note = "default CellChatDB; protein-mediated annotation categories were not found"
    ))
  }

  cellchat_db$interaction <- cellchat_db$interaction[keep, , drop = FALSE]
  return(list(database = cellchat_db, note = note))
}


SCassist_interactions_compute_centrality <- function(cellchat) {
  cellchat <- tryCatch(
    {
      CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    },
    error = function(e) {
      tryCatch(
        {
          CellChat::netAnalysis_computeCentrality(cellchat)
        },
        error = function(e2) {
          warning("CellChat centrality computation failed; centrality fields will be limited.")
          cellchat
        }
      )
    }
  )

  return(cellchat)
}


SCassist_interactions_build_metadata <- function(species,
                                                 database_info,
                                                 group_by,
                                                 assay,
                                                 average_method,
                                                 trim,
                                                 raw_use,
                                                 min_cells,
                                                 pval_threshold,
                                                 number_of_cells,
                                                 cell_group_sizes) {
  cell_group_size_values <- as.integer(cell_group_sizes)
  names(cell_group_size_values) <- names(cell_group_sizes)

  list(
    tool = "CellChat",
    cellchat_version = as.character(utils::packageVersion("CellChat")),
    species = species,
    database_used = database_info$database_name,
    database_scope = "default",
    database_note = database_info$database_note,
    group_by = group_by,
    assay = assay,
    average_method = average_method,
    trim = if (is.null(trim)) NA_real_ else trim,
    raw_use = raw_use,
    min_cells = min_cells,
    pval_threshold = pval_threshold,
    number_of_cells = number_of_cells,
    number_of_cell_groups = length(cell_group_sizes),
    cell_group_sizes = as.list(cell_group_size_values),
    cell_group_names = names(cell_group_sizes),
    generated_timestamp = as.character(Sys.time())
  )
}


SCassist_interactions_extract_communications <- function(cellchat,
                                                        cellchat_db,
                                                        pval_threshold) {
  communication <- tryCatch(
    {
      CellChat::subsetCommunication(cellchat)
    },
    error = function(e) {
      data.frame()
    }
  )

  if (is.null(communication) || nrow(communication) == 0) {
    return(SCassist_interactions_empty_interactions())
  }

  communication <- as.data.frame(communication, stringsAsFactors = FALSE)
  interactions <- SCassist_interactions_standardize_communication_columns(communication)
  interactions <- SCassist_interactions_attach_database_annotations(interactions, cellchat_db)
  interactions$communication_score <- suppressWarnings(as.numeric(interactions$communication_score))
  interactions$pval <- suppressWarnings(as.numeric(interactions$pval))

  if ("pval" %in% names(interactions) && any(!is.na(interactions$pval))) {
    interactions <- interactions[is.na(interactions$pval) | interactions$pval <= pval_threshold, , drop = FALSE]
  }

  if (nrow(interactions) == 0) {
    return(SCassist_interactions_empty_interactions())
  }

  interactions$autocrine <- as.character(interactions$source) == as.character(interactions$target)

  order_index <- order(interactions$communication_score, decreasing = TRUE, na.last = TRUE)
  interactions <- interactions[order_index, , drop = FALSE]
  interactions$rank_global <- seq_len(nrow(interactions))
  rownames(interactions) <- NULL

  return(interactions)
}


SCassist_interactions_empty_interactions <- function() {
  data.frame(
    source = character(),
    target = character(),
    ligand = character(),
    receptor = character(),
    interaction_name = character(),
    pathway_name = character(),
    communication_score = numeric(),
    pval = numeric(),
    annotation = character(),
    autocrine = logical(),
    rank_global = integer(),
    stringsAsFactors = FALSE
  )
}


SCassist_interactions_standardize_communication_columns <- function(communication) {
  source_col <- SCassist_interactions_find_column(communication, c("source", "sender", "from"))
  target_col <- SCassist_interactions_find_column(communication, c("target", "receiver", "to"))
  ligand_col <- SCassist_interactions_find_column(communication, c("ligand", "ligand.symbol", "ligand_symbol"))
  receptor_col <- SCassist_interactions_find_column(communication, c("receptor", "receptor.symbol", "receptor_symbol"))
  interaction_col <- SCassist_interactions_find_column(communication, c("interaction_name", "interaction_name_2", "interaction"))
  pathway_col <- SCassist_interactions_find_column(communication, c("pathway_name", "pathway", "signaling"))
  prob_col <- SCassist_interactions_find_column(communication, c("prob", "communication_score", "score", "weight"))
  pval_col <- SCassist_interactions_find_column(communication, c("pval", "p.value", "p_value", "p"))
  annotation_col <- SCassist_interactions_find_column(communication, c("annotation", "category", "interaction_category"))

  data.frame(
    source = SCassist_interactions_get_column_or_na(communication, source_col),
    target = SCassist_interactions_get_column_or_na(communication, target_col),
    ligand = SCassist_interactions_get_column_or_na(communication, ligand_col),
    receptor = SCassist_interactions_get_column_or_na(communication, receptor_col),
    interaction_name = SCassist_interactions_get_column_or_na(communication, interaction_col),
    pathway_name = SCassist_interactions_get_column_or_na(communication, pathway_col),
    communication_score = SCassist_interactions_get_column_or_na(communication, prob_col),
    pval = SCassist_interactions_get_column_or_na(communication, pval_col),
    annotation = SCassist_interactions_get_column_or_na(communication, annotation_col),
    stringsAsFactors = FALSE
  )
}


SCassist_interactions_find_column <- function(data, candidates) {
  matches <- candidates[candidates %in% names(data)]
  if (length(matches) > 0) {
    return(matches[[1]])
  }
  return(NULL)
}


SCassist_interactions_get_column_or_na <- function(data, column_name) {
  if (is.null(column_name)) {
    return(rep(NA, nrow(data)))
  }
  return(data[[column_name]])
}


SCassist_interactions_attach_database_annotations <- function(interactions, cellchat_db) {
  if (!"annotation" %in% names(interactions)) {
    interactions$annotation <- NA_character_
  }

  has_annotation <- any(!is.na(interactions$annotation) & interactions$annotation != "")
  if (has_annotation || is.null(cellchat_db$interaction)) {
    return(interactions)
  }

  db_interaction <- as.data.frame(cellchat_db$interaction, stringsAsFactors = FALSE)
  db_interaction_col <- SCassist_interactions_find_column(
    db_interaction,
    c("interaction_name", "interaction_name_2", "interaction")
  )
  db_annotation_col <- SCassist_interactions_find_column(
    db_interaction,
    c("annotation", "category", "interaction_category")
  )

  if (is.null(db_interaction_col) || is.null(db_annotation_col) ||
      !"interaction_name" %in% names(interactions)) {
    return(interactions)
  }

  match_index <- match(interactions$interaction_name, db_interaction[[db_interaction_col]])
  interactions$annotation <- db_interaction[[db_annotation_col]][match_index]
  return(interactions)
}


SCassist_interactions_summarize_pathways <- function(interactions, cell_groups = NULL) {
  if (nrow(interactions) == 0 || !"pathway_name" %in% names(interactions)) {
    return(SCassist_interactions_empty_pathway_summary())
  }

  pathway_names <- unique(interactions$pathway_name)
  pathway_names <- pathway_names[!is.na(pathway_names) & pathway_names != ""]

  if (length(pathway_names) == 0) {
    return(SCassist_interactions_empty_pathway_summary())
  }

  pathway_summary <- lapply(pathway_names, function(pathway_name) {
    pathway_rows <- interactions[interactions$pathway_name == pathway_name, , drop = FALSE]
    score <- suppressWarnings(as.numeric(pathway_rows$communication_score))
    score[is.na(score)] <- 0
    total_score <- sum(score)
    sender_scores <- SCassist_interactions_sum_by(pathway_rows$source, score)
    receiver_scores <- SCassist_interactions_sum_by(pathway_rows$target, score)
    top_lr <- SCassist_interactions_top_lr(pathway_rows, score, n = 5)
    autocrine_count <- sum(pathway_rows$autocrine, na.rm = TRUE)
    paracrine_count <- nrow(pathway_rows) - autocrine_count
    topology_label <- SCassist_interactions_topology_label(sender_scores, receiver_scores)

    data.frame(
      pathway_name = pathway_name,
      number_of_significant_interactions = nrow(pathway_rows),
      total_communication_score = total_score,
      mean_communication_score = if (length(score) > 0) mean(score) else NA_real_,
      dominant_senders = SCassist_interactions_collapse_top_names(sender_scores, n = 3),
      dominant_receivers = SCassist_interactions_collapse_top_names(receiver_scores, n = 3),
      top_ligand_receptor_interactions = paste(top_lr, collapse = "; "),
      autocrine_count = autocrine_count,
      paracrine_count = paracrine_count,
      autocrine_fraction = if (nrow(pathway_rows) > 0) autocrine_count / nrow(pathway_rows) else NA_real_,
      paracrine_fraction = if (nrow(pathway_rows) > 0) paracrine_count / nrow(pathway_rows) else NA_real_,
      topology_label = topology_label,
      stringsAsFactors = FALSE
    )
  })

  pathway_summary <- do.call(rbind, pathway_summary)
  pathway_summary <- pathway_summary[order(pathway_summary$total_communication_score, decreasing = TRUE), , drop = FALSE]
  rownames(pathway_summary) <- NULL

  return(pathway_summary)
}


SCassist_interactions_empty_pathway_summary <- function() {
  data.frame(
    pathway_name = character(),
    number_of_significant_interactions = integer(),
    total_communication_score = numeric(),
    mean_communication_score = numeric(),
    dominant_senders = character(),
    dominant_receivers = character(),
    top_ligand_receptor_interactions = character(),
    autocrine_count = integer(),
    paracrine_count = integer(),
    autocrine_fraction = numeric(),
    paracrine_fraction = numeric(),
    topology_label = character(),
    stringsAsFactors = FALSE
  )
}


SCassist_interactions_summarize_roles <- function(interactions, cell_groups, cellchat = NULL) {
  if (length(cell_groups) == 0) {
    cell_groups <- sort(unique(c(interactions$source, interactions$target)))
  }

  outgoing <- SCassist_interactions_score_by_group(interactions, "source", cell_groups)
  incoming <- SCassist_interactions_score_by_group(interactions, "target", cell_groups)
  centrality <- SCassist_interactions_extract_centrality(cellchat, cell_groups)

  role_summary <- data.frame(
    cell_group = cell_groups,
    outgoing_strength = as.numeric(outgoing[cell_groups]),
    incoming_strength = as.numeric(incoming[cell_groups]),
    outdegree_rank = rank(-as.numeric(outgoing[cell_groups]), ties.method = "min"),
    indegree_rank = rank(-as.numeric(incoming[cell_groups]), ties.method = "min"),
    flow_betweenness = centrality$flow_betweenness,
    information_centrality = centrality$information_centrality,
    top_outgoing_pathways = vapply(cell_groups, function(group_name) {
      SCassist_interactions_top_pathways_for_group(interactions, group_name, "source")
    }, character(1)),
    top_incoming_pathways = vapply(cell_groups, function(group_name) {
      SCassist_interactions_top_pathways_for_group(interactions, group_name, "target")
    }, character(1)),
    stringsAsFactors = FALSE
  )

  role_summary$role_label <- SCassist_interactions_assign_role_labels(role_summary)
  role_summary <- role_summary[order(role_summary$outdegree_rank, role_summary$indegree_rank), , drop = FALSE]
  rownames(role_summary) <- NULL

  return(role_summary)
}


SCassist_interactions_score_by_group <- function(interactions, column_name, cell_groups) {
  scores <- stats::setNames(rep(0, length(cell_groups)), cell_groups)

  if (nrow(interactions) == 0 || !column_name %in% names(interactions)) {
    return(scores)
  }

  score <- suppressWarnings(as.numeric(interactions$communication_score))
  score[is.na(score)] <- 0
  group_scores <- SCassist_interactions_sum_by(interactions[[column_name]], score)
  scores[names(group_scores)] <- group_scores

  return(scores)
}


SCassist_interactions_extract_centrality <- function(cellchat, cell_groups) {
  result <- list(
    flow_betweenness = rep(NA_real_, length(cell_groups)),
    information_centrality = rep(NA_real_, length(cell_groups))
  )

  if (is.null(cellchat)) {
    return(result)
  }

  centr <- tryCatch(
    {
      cellchat@netP$centr
    },
    error = function(e) {
      NULL
    }
  )

  if (is.null(centr) || length(centr) == 0) {
    return(result)
  }

  flow_values <- stats::setNames(rep(0, length(cell_groups)), cell_groups)
  info_values <- stats::setNames(rep(0, length(cell_groups)), cell_groups)
  flow_seen <- FALSE
  info_seen <- FALSE

  for (centr_item in centr) {
    extracted <- SCassist_interactions_extract_centrality_item(centr_item, cell_groups)

    if (!all(is.na(extracted$flow_betweenness))) {
      flow_part <- extracted$flow_betweenness
      flow_part[is.na(flow_part)] <- 0
      flow_values <- flow_values + flow_part
      flow_seen <- TRUE
    }

    if (!all(is.na(extracted$information_centrality))) {
      info_part <- extracted$information_centrality
      info_part[is.na(info_part)] <- 0
      info_values <- info_values + info_part
      info_seen <- TRUE
    }
  }

  if (flow_seen) {
    result$flow_betweenness <- as.numeric(flow_values[cell_groups])
  }

  if (info_seen) {
    result$information_centrality <- as.numeric(info_values[cell_groups])
  }

  return(result)
}


SCassist_interactions_extract_centrality_item <- function(centr_item, cell_groups) {
  flow <- stats::setNames(rep(NA_real_, length(cell_groups)), cell_groups)
  info <- stats::setNames(rep(NA_real_, length(cell_groups)), cell_groups)

  if (is.data.frame(centr_item) || is.matrix(centr_item)) {
    centr_df <- as.data.frame(centr_item)
    row_names <- rownames(centr_df)
    flow_col <- SCassist_interactions_find_column(
      centr_df,
      c("flowbet", "flow_betweenness", "betweenness", "flow.betweenness")
    )
    info_col <- SCassist_interactions_find_column(
      centr_df,
      c("info", "information", "information_centrality", "information.centrality")
    )

    if (!is.null(flow_col)) {
      values <- suppressWarnings(as.numeric(centr_df[[flow_col]]))
      names(values) <- row_names
      flow[names(values)] <- values
    }

    if (!is.null(info_col)) {
      values <- suppressWarnings(as.numeric(centr_df[[info_col]]))
      names(values) <- row_names
      info[names(values)] <- values
    }
  } else if (is.list(centr_item)) {
    flow_candidates <- c("flowbet", "flow_betweenness", "betweenness", "flow.betweenness")
    info_candidates <- c("info", "information", "information_centrality", "information.centrality")

    flow_name <- flow_candidates[flow_candidates %in% names(centr_item)][1]
    info_name <- info_candidates[info_candidates %in% names(centr_item)][1]

    if (!is.na(flow_name) && !is.null(centr_item[[flow_name]])) {
      values <- suppressWarnings(as.numeric(centr_item[[flow_name]]))
      names(values) <- names(centr_item[[flow_name]])
      flow[names(values)] <- values
    }

    if (!is.na(info_name) && !is.null(centr_item[[info_name]])) {
      values <- suppressWarnings(as.numeric(centr_item[[info_name]]))
      names(values) <- names(centr_item[[info_name]])
      info[names(values)] <- values
    }
  }

  list(
    flow_betweenness = as.numeric(flow[cell_groups]),
    information_centrality = as.numeric(info[cell_groups])
  )
}


SCassist_interactions_assign_role_labels <- function(role_summary) {
  outgoing <- role_summary$outgoing_strength
  incoming <- role_summary$incoming_strength
  total_strength <- outgoing + incoming

  if (length(total_strength) == 0 || all(is.na(total_strength) | total_strength == 0)) {
    return(rep("low_connectivity", nrow(role_summary)))
  }

  outgoing_high <- outgoing >= stats::quantile(outgoing, probs = 0.75, na.rm = TRUE)
  incoming_high <- incoming >= stats::quantile(incoming, probs = 0.75, na.rm = TRUE)
  total_mid <- total_strength >= stats::median(total_strength, na.rm = TRUE)

  labels <- rep("low_connectivity", nrow(role_summary))
  labels[outgoing_high & !incoming_high] <- "dominant_sender"
  labels[incoming_high & !outgoing_high] <- "dominant_receiver"
  labels[outgoing_high & incoming_high] <- "bidirectional_hub"

  flow <- role_summary$flow_betweenness
  flow_high <- !is.na(flow) & flow >= stats::quantile(flow, probs = 0.75, na.rm = TRUE)
  labels[flow_high & !(labels %in% c("dominant_sender", "dominant_receiver", "bidirectional_hub"))] <- "mediator"
  labels[total_mid & labels == "low_connectivity"] <- "influencer"

  return(labels)
}


SCassist_interactions_summarize_autocrine_paracrine <- function(interactions) {
  if (nrow(interactions) == 0) {
    return(data.frame(
      total_autocrine_interactions = 0,
      total_paracrine_interactions = 0,
      autocrine_probability_mass = 0,
      paracrine_probability_mass = 0,
      top_autocrine_pathways = NA_character_,
      top_paracrine_pathways = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  score <- suppressWarnings(as.numeric(interactions$communication_score))
  score[is.na(score)] <- 0

  autocrine_rows <- interactions$autocrine
  paracrine_rows <- !interactions$autocrine

  data.frame(
    total_autocrine_interactions = sum(autocrine_rows, na.rm = TRUE),
    total_paracrine_interactions = sum(paracrine_rows, na.rm = TRUE),
    autocrine_probability_mass = sum(score[autocrine_rows], na.rm = TRUE),
    paracrine_probability_mass = sum(score[paracrine_rows], na.rm = TRUE),
    top_autocrine_pathways = SCassist_interactions_top_pathways(interactions[autocrine_rows, , drop = FALSE], n = 5),
    top_paracrine_pathways = SCassist_interactions_top_pathways(interactions[paracrine_rows, , drop = FALSE], n = 5),
    stringsAsFactors = FALSE
  )
}


SCassist_interactions_summarize_interaction_types <- function(interactions) {
  if (nrow(interactions) == 0 || !"annotation" %in% names(interactions) ||
      all(is.na(interactions$annotation) | interactions$annotation == "")) {
    return(data.frame(
      annotation = NA_character_,
      number_of_interactions = 0,
      total_communication_score = NA_real_,
      mean_communication_score = NA_real_,
      note = "CellChatDB interaction annotations were unavailable in the extracted communication table.",
      stringsAsFactors = FALSE
    ))
  }

  annotations <- unique(interactions$annotation)
  annotations <- annotations[!is.na(annotations) & annotations != ""]

  type_summary <- lapply(annotations, function(annotation) {
    rows <- interactions[interactions$annotation == annotation, , drop = FALSE]
    score <- suppressWarnings(as.numeric(rows$communication_score))
    score[is.na(score)] <- 0
    data.frame(
      annotation = annotation,
      number_of_interactions = nrow(rows),
      total_communication_score = sum(score),
      mean_communication_score = if (length(score) > 0) mean(score) else NA_real_,
      note = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  type_summary <- do.call(rbind, type_summary)
  type_summary <- type_summary[order(type_summary$total_communication_score, decreasing = TRUE), , drop = FALSE]
  rownames(type_summary) <- NULL

  return(type_summary)
}


SCassist_interactions_build_llm_context <- function(metadata,
                                                    interactions,
                                                    pathway_summary,
                                                    cell_role_summary,
                                                    autocrine_paracrine_summary,
                                                    interaction_type_summary,
                                                    top_n_interactions,
                                                    top_n_pathways,
                                                    experimental_context = NULL) {
  top_interactions <- head(interactions, top_n_interactions)
  top_pathways <- head(pathway_summary, top_n_pathways)

  list(
    metadata = metadata,
    global_network_summary = SCassist_interactions_global_summary(
      interactions = interactions,
      pathway_summary = pathway_summary,
      cell_role_summary = cell_role_summary,
      autocrine_paracrine_summary = autocrine_paracrine_summary
    ),
    top_pathways = top_pathways,
    top_interactions = top_interactions,
    cell_role_summary = cell_role_summary,
    autocrine_paracrine_summary = autocrine_paracrine_summary,
    interaction_type_summary = interaction_type_summary,
    experimental_context = experimental_context,
    interpretation_caveats = c(
      "CellChat infers potential communication from mRNA expression and curated prior knowledge.",
      "CellChat communication probability is an interaction strength score, not a literal probability.",
      "CellChat does not prove physical signaling or causality.",
      "Protein-level activity and spatial proximity are not directly measured from standard scRNA-seq.",
      "Cell group definitions strongly affect results.",
      "triMean is conservative and tends to produce fewer but stronger interactions.",
      "Weak or missing signaling may depend on sequencing depth, average expression method, and cell group size."
    )
  )
}


SCassist_interactions_global_summary <- function(interactions,
                                                pathway_summary,
                                                cell_role_summary,
                                                autocrine_paracrine_summary) {
  strongest_sender <- if (nrow(cell_role_summary) > 0) {
    cell_role_summary$cell_group[which.max(cell_role_summary$outgoing_strength)]
  } else {
    NA_character_
  }

  strongest_receiver <- if (nrow(cell_role_summary) > 0) {
    cell_role_summary$cell_group[which.max(cell_role_summary$incoming_strength)]
  } else {
    NA_character_
  }

  dominant_pathway <- if (nrow(pathway_summary) > 0) {
    pathway_summary$pathway_name[which.max(pathway_summary$total_communication_score)]
  } else {
    NA_character_
  }

  list(
    total_significant_interactions = nrow(interactions),
    number_of_significant_pathways = nrow(pathway_summary),
    total_communication_score = if (nrow(interactions) > 0) {
      sum(suppressWarnings(as.numeric(interactions$communication_score)), na.rm = TRUE)
    } else {
      0
    },
    strongest_sender = strongest_sender,
    strongest_receiver = strongest_receiver,
    dominant_pathway = dominant_pathway,
    autocrine_interactions = autocrine_paracrine_summary$total_autocrine_interactions,
    paracrine_interactions = autocrine_paracrine_summary$total_paracrine_interactions,
    autocrine_probability_mass = autocrine_paracrine_summary$autocrine_probability_mass,
    paracrine_probability_mass = autocrine_paracrine_summary$paracrine_probability_mass
  )
}


SCassist_interactions_build_llm_prompt <- function(llm_context) {
  context_text <- SCassist_interactions_context_to_text(llm_context)

  paste0(
    "You are SCassist InteractionAgent, an AI assistant for interpreting CellChat v2 single-condition cell-cell communication analysis.\n\n",
    "Interpret the CellChat results biologically using cautious scientific wording. ",
    "Use phrases such as 'CellChat inferred', 'The data suggest', and 'This is consistent with'. ",
    "Avoid phrases such as 'This proves' or 'These cells definitely communicate'.\n\n",
    "Your response should include:\n",
    "1. Overall communication landscape\n",
    "2. Dominant sender cell groups\n",
    "3. Dominant receiver cell groups\n",
    "4. Major signaling pathways\n",
    "5. Top ligand-receptor axes\n",
    "6. Autocrine versus paracrine communication\n",
    "7. Biological interpretation of likely processes\n",
    "8. Confidence/caveats\n",
    "9. Suggested follow-up validation or analyses\n\n",
    "Use the structured context below. Do not treat the communication score as a literal probability, ",
    "and do not infer causality or physical contact without validation.\n\n",
    "Structured CellChat context:\n",
    context_text
  )
}


SCassist_interactions_context_to_text <- function(llm_context) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(paste0(
      "```json\n",
      jsonlite::toJSON(
        llm_context,
        pretty = TRUE,
        auto_unbox = TRUE,
        na = "null"
      ),
      "\n```"
    ))
  }

  paste(capture.output(str(llm_context, max.level = 4)), collapse = "\n")
}


SCassist_interactions_add_llm_response <- function(interaction_results,
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
    interaction_results["llm_response"] <- list(NULL)
    return(interaction_results)
  }

  if (llm_server == "google") {
    interaction_results$llm_response <- SCassist_analyze_interactions_G(
      llm_prompt = interaction_results$llm_prompt,
      temperature = temperature,
      max_output_tokens = max_output_tokens,
      model = model_G,
      api_key_file = api_key_file
    )
  } else if (llm_server == "ollama") {
    interaction_results$llm_response <- SCassist_analyze_interactions_L(
      llm_prompt = interaction_results$llm_prompt,
      model_params = model_params,
      model = model_O
    )
  } else if (llm_server == "openai") {
    interaction_results$llm_response <- SCassist_analyze_interactions_C(
      llm_prompt = interaction_results$llm_prompt,
      temperature = temperature,
      max_output_tokens = max_output_tokens,
      model = model_C,
      api_key_file = api_key_file
    )
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama' or 'openai'.", call. = FALSE)
  }

  return(interaction_results)
}


SCassist_analyze_interactions_G <- function(llm_prompt,
                                            temperature = 0,
                                            max_output_tokens = 10048,
                                            model = "gemini-1.5-flash-latest",
                                            api_key_file = "api_keys.txt") {
  model_query <- paste0(model, ":generateContent")
  api_key <- readLines(api_key_file)

  response <- httr::POST(
    url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
    query = list(key = api_key),
    httr::content_type_json(),
    encode = "json",
    body = list(
      contents = list(
        parts = list(
          list(text = llm_prompt)
        )
      ),
      generationConfig = list(
        temperature = temperature,
        maxOutputTokens = max_output_tokens,
        seed = 123456
      )
    )
  )

  if (response$status_code > 200) {
    stop(paste("Error - ", httr::content(response)$error$message))
  }

  candidates <- httr::content(response)$candidates
  outputs <- unlist(lapply(candidates, function(candidate) candidate$content$parts))
  generated_text <- outputs[["text"]]

  cat(generated_text)
  return(generated_text)
}


SCassist_analyze_interactions_L <- function(llm_prompt,
                                            model_params = list(seed = 42, temperature = 0, num_gpu = 0),
                                            model = "llama3") {
  response <- tryCatch(
    {
      suppressMessages(
        rollama::query(
          llm_prompt,
          model = model,
          screen = FALSE,
          model_params = model_params
        )
      )
    },
    error = function(e) {
      stop("Error: The LLM query encountered an error. Please check your ollama server connection and/or the model.", call. = FALSE)
    }
  )

  if (is.null(response[[1]]$message$content)) {
    stop("Error: The LLM returned an invalid response. Please check the LLM model and parameters.", call. = FALSE)
  }

  cat(response[[1]]$message$content)
  return(response[[1]]$message$content)
}


SCassist_analyze_interactions_C <- function(llm_prompt,
                                            temperature = 0,
                                            max_output_tokens = 10048,
                                            model = "gpt-4o-mini",
                                            api_key_file = "api_keys.txt") {
  api_key <- readLines(api_key_file, encoding = "UTF-8")

  response <- httr::POST(
    url = "https://api.openai.com/v1/chat/completions",
    httr::add_headers("Authorization" = paste("Bearer", api_key)),
    httr::content_type_json(),
    encode = "json",
    body = list(
      model = model,
      messages = list(
        list(
          role = "user",
          content = llm_prompt
        )
      ),
      temperature = temperature,
      max_tokens = max_output_tokens
    )
  )

  if (httr::http_error(response)) {
    stop(paste(
      "OpenAI API request failed with status",
      httr::http_status(response)$message,
      "\nContent:",
      httr::content(response, "text")
    ))
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required but not installed.", call. = FALSE)
  }

  response_content <- httr::content(response, "text")
  response_json <- jsonlite::fromJSON(response_content, flatten = TRUE)

  if (!is.null(response_json$choices) &&
      length(response_json$choices) > 0 &&
      !is.null(response_json$choices$message.content)) {
    generated_text <- response_json$choices$message.content
    cat(generated_text)
  } else {
    cat("Error: Could not extract generated text from OpenAI response.\n")
    print(response_json)
    generated_text <- NULL
  }

  return(generated_text)
}


SCassist_interactions_write_outputs <- function(interaction_results,
                                                output_dir,
                                                run_llm) {
  if (is.null(output_dir)) {
    return(invisible(NULL))
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  saveRDS(interaction_results$cellchat_object, file = file.path(output_dir, "cellchat_object.rds"))
  utils::write.csv(interaction_results$interactions, file = file.path(output_dir, "cellchat_interactions.csv"), row.names = FALSE)
  utils::write.csv(interaction_results$pathway_summary, file = file.path(output_dir, "cellchat_pathway_summary.csv"), row.names = FALSE)
  utils::write.csv(interaction_results$cell_role_summary, file = file.path(output_dir, "cellchat_cell_role_summary.csv"), row.names = FALSE)
  utils::write.csv(interaction_results$autocrine_paracrine_summary, file = file.path(output_dir, "cellchat_autocrine_paracrine_summary.csv"), row.names = FALSE)
  utils::write.csv(interaction_results$interaction_type_summary, file = file.path(output_dir, "cellchat_interaction_type_summary.csv"), row.names = FALSE)

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(
      interaction_results$llm_context,
      path = file.path(output_dir, "cellchat_llm_context.json"),
      pretty = TRUE,
      auto_unbox = TRUE,
      na = "null"
    )
  } else {
    warning("Package 'jsonlite' is not available; skipping cellchat_llm_context.json output.")
  }

  if (run_llm && !is.null(interaction_results$llm_response)) {
    writeLines(interaction_results$llm_response, con = file.path(output_dir, "cellchat_interpretation.txt"))
  }

  return(invisible(NULL))
}


SCassist_interactions_sum_by <- function(groups, values) {
  groups <- as.character(groups)
  valid <- !is.na(groups) & groups != ""
  groups <- groups[valid]
  values <- values[valid]

  if (length(groups) == 0) {
    return(numeric())
  }

  result <- tapply(values, groups, sum, na.rm = TRUE)
  result <- sort(result, decreasing = TRUE)
  return(result)
}


SCassist_interactions_collapse_top_names <- function(named_scores, n = 3) {
  if (length(named_scores) == 0) {
    return(NA_character_)
  }

  top_names <- names(head(sort(named_scores, decreasing = TRUE), n))
  paste(top_names, collapse = ", ")
}


SCassist_interactions_top_lr <- function(pathway_rows, score, n = 5) {
  if (nrow(pathway_rows) == 0) {
    return(character())
  }

  lr <- if ("ligand" %in% names(pathway_rows) && "receptor" %in% names(pathway_rows)) {
    paste(pathway_rows$ligand, pathway_rows$receptor, sep = "-")
  } else {
    pathway_rows$interaction_name
  }

  lr[is.na(lr) | lr == "NA-NA" | lr == ""] <- pathway_rows$interaction_name[is.na(lr) | lr == "NA-NA" | lr == ""]
  order_index <- order(score, decreasing = TRUE, na.last = TRUE)
  unique(lr[order_index])[seq_len(min(n, length(unique(lr[order_index]))))]
}


SCassist_interactions_topology_label <- function(sender_scores, receiver_scores) {
  sender_focused <- SCassist_interactions_is_focused(sender_scores)
  receiver_focused <- SCassist_interactions_is_focused(receiver_scores)

  if (sender_focused && receiver_focused) {
    return("focused_sender_focused_receiver")
  }
  if (sender_focused && !receiver_focused) {
    return("focused_sender_broad_receiver")
  }
  if (!sender_focused && receiver_focused) {
    return("broad_sender_focused_receiver")
  }
  return("broad_sender_broad_receiver")
}


SCassist_interactions_is_focused <- function(scores) {
  if (length(scores) <= 1) {
    return(TRUE)
  }

  total <- sum(scores, na.rm = TRUE)
  if (total == 0) {
    return(FALSE)
  }

  max(scores, na.rm = TRUE) / total >= 0.5
}


SCassist_interactions_top_pathways <- function(interactions, n = 5) {
  if (nrow(interactions) == 0 || !"pathway_name" %in% names(interactions)) {
    return(NA_character_)
  }

  score <- suppressWarnings(as.numeric(interactions$communication_score))
  score[is.na(score)] <- 0
  pathway_scores <- SCassist_interactions_sum_by(interactions$pathway_name, score)

  if (length(pathway_scores) == 0) {
    return(NA_character_)
  }

  paste(names(head(pathway_scores, n)), collapse = ", ")
}


SCassist_interactions_top_pathways_for_group <- function(interactions, group_name, column_name, n = 5) {
  if (nrow(interactions) == 0 || !column_name %in% names(interactions)) {
    return(NA_character_)
  }

  rows <- interactions[as.character(interactions[[column_name]]) == as.character(group_name), , drop = FALSE]
  SCassist_interactions_top_pathways(rows, n = n)
}
