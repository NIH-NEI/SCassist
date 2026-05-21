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


