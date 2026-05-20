library(testthat)
library(SCassist)

test_that("SCassist_analyze_interactions handles invalid object names", {
  expect_error(
    SCassist_analyze_interactions("nonexistent_object", run_llm = FALSE),
    "not found"
  )
})

test_that("SCassist_analyze_interactions validates species", {
  expect_error(
    SCassist_analyze_interactions("nonexistent_object", species = "rat", run_llm = FALSE),
    "species"
  )
})

test_that("SCassist_analyze_interactions validates group_by metadata", {
  skip_if_not_installed("Seurat")

  counts <- matrix(
    c(1, 0, 2, 0, 3, 1, 0, 2, 1, 4, 0, 1),
    nrow = 3,
    dimnames = list(c("GeneA", "GeneB", "GeneC"), paste0("cell", 1:4))
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)
  assign("scassist_interactions_missing_group", obj, envir = .GlobalEnv)
  on.exit(rm("scassist_interactions_missing_group", envir = .GlobalEnv), add = TRUE)

  expect_error(
    SCassist_analyze_interactions(
      "scassist_interactions_missing_group",
      group_by = "celltype",
      run_llm = FALSE
    ),
    "group_by"
  )
})

test_that("run_llm FALSE does not require an API key", {
  interaction_results <- list(
    llm_prompt = "Interpret this compact CellChat context.",
    llm_response = "previous response"
  )

  result <- SCassist:::SCassist_interactions_add_llm_response(
    interaction_results = interaction_results,
    run_llm = FALSE,
    llm_server = "google",
    temperature = 0,
    max_output_tokens = 10048,
    model_G = "gemini-1.5-flash-latest",
    model_O = "llama3",
    model_C = "gpt-4o-mini",
    api_key_file = "missing_api_key_file.txt",
    model_params = list()
  )

  expect_null(result$llm_response)
  expect_true("llm_response" %in% names(result))
})

test_that("interaction summary helpers work on small synthetic tables", {
  interactions <- data.frame(
    source = c("T cell", "T cell", "B cell", "B cell"),
    target = c("B cell", "T cell", "B cell", "T cell"),
    ligand = c("L1", "L2", "L3", "L4"),
    receptor = c("R1", "R2", "R3", "R4"),
    interaction_name = c("L1_R1", "L2_R2", "L3_R3", "L4_R4"),
    pathway_name = c("CXCL", "CXCL", "MIF", "MIF"),
    communication_score = c(0.8, 0.4, 0.6, 0.2),
    pval = c(0.01, 0.02, 0.03, 0.04),
    annotation = c("Secreted Signaling", "Secreted Signaling", "ECM-Receptor", "ECM-Receptor"),
    autocrine = c(FALSE, TRUE, TRUE, FALSE),
    rank_global = 1:4,
    stringsAsFactors = FALSE
  )

  pathway_summary <- SCassist:::SCassist_interactions_summarize_pathways(
    interactions,
    cell_groups = c("T cell", "B cell")
  )
  cell_role_summary <- SCassist:::SCassist_interactions_summarize_roles(
    interactions,
    cell_groups = c("T cell", "B cell"),
    cellchat = NULL
  )
  ap_summary <- SCassist:::SCassist_interactions_summarize_autocrine_paracrine(interactions)
  type_summary <- SCassist:::SCassist_interactions_summarize_interaction_types(interactions)

  expect_equal(nrow(pathway_summary), 2)
  expect_true("top_ligand_receptor_interactions" %in% names(pathway_summary))
  expect_equal(nrow(cell_role_summary), 2)
  expect_true("role_label" %in% names(cell_role_summary))
  expect_equal(ap_summary$total_autocrine_interactions, 2)
  expect_equal(ap_summary$total_paracrine_interactions, 2)
  expect_equal(nrow(type_summary), 2)
})

test_that("LLM context and prompt contain CellChat caveats", {
  interactions <- data.frame(
    source = "T cell",
    target = "B cell",
    ligand = "L1",
    receptor = "R1",
    interaction_name = "L1_R1",
    pathway_name = "CXCL",
    communication_score = 0.8,
    pval = 0.01,
    annotation = "Secreted Signaling",
    autocrine = FALSE,
    rank_global = 1,
    stringsAsFactors = FALSE
  )

  metadata <- list(
    tool = "CellChat",
    species = "human",
    group_by = "celltype",
    average_method = "triMean"
  )
  pathway_summary <- SCassist:::SCassist_interactions_summarize_pathways(interactions, c("T cell", "B cell"))
  cell_role_summary <- SCassist:::SCassist_interactions_summarize_roles(interactions, c("T cell", "B cell"), NULL)
  ap_summary <- SCassist:::SCassist_interactions_summarize_autocrine_paracrine(interactions)
  type_summary <- SCassist:::SCassist_interactions_summarize_interaction_types(interactions)

  llm_context <- SCassist:::SCassist_interactions_build_llm_context(
    metadata = metadata,
    interactions = interactions,
    pathway_summary = pathway_summary,
    cell_role_summary = cell_role_summary,
    autocrine_paracrine_summary = ap_summary,
    interaction_type_summary = type_summary,
    top_n_interactions = 10,
    top_n_pathways = 10,
    experimental_context = "synthetic test"
  )
  llm_prompt <- SCassist:::SCassist_interactions_build_llm_prompt(llm_context)

  expect_true("interpretation_caveats" %in% names(llm_context))
  expect_match(llm_prompt, "CellChat inferred")
  expect_match(llm_prompt, "not a literal probability")
})

test_that("CellChat database helper loads installed package data", {
  skip_if_not_installed("CellChat")

  db_info <- SCassist:::SCassist_interactions_get_cellchat_db(
    species = "human",
    database_scope = "default"
  )

  expect_equal(db_info$database_name, "CellChatDB.human")
  expect_true("database_note" %in% names(db_info))
  expect_true(is.list(db_info$database))
  expect_true("interaction" %in% names(db_info$database))
  expect_gt(nrow(db_info$database$interaction), 0)
})

test_that("missing CellChat dependency reports an installation hint", {
  skip_if_not_installed("Seurat")
  if (requireNamespace("CellChat", quietly = TRUE)) {
    skip("CellChat is installed; missing-dependency branch is not applicable.")
  }

  counts <- matrix(
    c(1, 0, 2, 0, 3, 1, 0, 2, 1, 4, 0, 1),
    nrow = 3,
    dimnames = list(c("GeneA", "GeneB", "GeneC"), paste0("cell", 1:4))
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)
  obj$celltype <- c("A", "A", "B", "B")
  assign("scassist_interactions_cellchat_missing", obj, envir = .GlobalEnv)
  on.exit(rm("scassist_interactions_cellchat_missing", envir = .GlobalEnv), add = TRUE)

  expect_error(
    SCassist_analyze_interactions(
      "scassist_interactions_cellchat_missing",
      run_llm = FALSE
    ),
    "jinworks/CellChat"
  )
})
