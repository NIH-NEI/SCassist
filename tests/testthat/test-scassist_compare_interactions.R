library(testthat)
library(SCassist)

scassist_compare_test_object <- function(groups = c("Sender", "Sender", "Receiver", "Receiver", "Sender", "Sender", "Receiver", "Receiver"),
                                         conditions = c("control", "control", "control", "control", "treated", "treated", "treated", "treated")) {
  skip_if_not_installed("Seurat")

  counts <- matrix(
    c(
      5, 1, 0, 0, 4, 1, 0, 0,
      0, 4, 1, 0, 0, 5, 1, 0,
      1, 0, 5, 1, 1, 0, 6, 1,
      0, 1, 0, 4, 0, 1, 0, 5,
      2, 2, 2, 2, 3, 3, 3, 3
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:5), paste0("cell", seq_along(groups)))
  )

  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  obj$celltype <- groups
  obj$condition <- conditions
  obj
}

test_that("SCassist_compare_interactions handles invalid object names", {
  expect_error(
    SCassist_compare_interactions(
      "nonexistent_object",
      condition_a = "control",
      condition_b = "treated",
      run_llm = FALSE
    ),
    "not found"
  )
})

test_that("SCassist_compare_interactions validates species", {
  expect_error(
    SCassist_compare_interactions(
      "nonexistent_object",
      condition_a = "control",
      condition_b = "treated",
      species = "rat",
      run_llm = FALSE
    ),
    "species"
  )
})

test_that("SCassist_compare_interactions validates missing condition arguments", {
  expect_error(
    SCassist_compare_interactions(
      "nonexistent_object",
      condition_b = "treated",
      run_llm = FALSE
    ),
    "condition_a"
  )

  expect_error(
    SCassist_compare_interactions(
      "nonexistent_object",
      condition_a = "control",
      run_llm = FALSE
    ),
    "condition_b"
  )
})

test_that("SCassist_compare_interactions validates distinct conditions", {
  expect_error(
    SCassist_compare_interactions(
      "nonexistent_object",
      condition_a = "control",
      condition_b = "control",
      run_llm = FALSE
    ),
    "different"
  )
})

test_that("SCassist_compare_interactions validates group_by metadata", {
  obj <- scassist_compare_test_object()
  obj$celltype <- NULL
  assign("scassist_compare_missing_group", obj, envir = .GlobalEnv)
  on.exit(rm("scassist_compare_missing_group", envir = .GlobalEnv), add = TRUE)

  expect_error(
    SCassist_compare_interactions(
      "scassist_compare_missing_group",
      group_by = "celltype",
      condition_by = "condition",
      condition_a = "control",
      condition_b = "treated",
      run_llm = FALSE
    ),
    "group_by"
  )
})

test_that("SCassist_compare_interactions validates condition_by metadata", {
  obj <- scassist_compare_test_object()
  obj$condition <- NULL
  assign("scassist_compare_missing_condition", obj, envir = .GlobalEnv)
  on.exit(rm("scassist_compare_missing_condition", envir = .GlobalEnv), add = TRUE)

  expect_error(
    SCassist_compare_interactions(
      "scassist_compare_missing_condition",
      group_by = "celltype",
      condition_by = "condition",
      condition_a = "control",
      condition_b = "treated",
      run_llm = FALSE
    ),
    "condition_by"
  )
})

test_that("SCassist_compare_interactions validates requested condition values", {
  obj <- scassist_compare_test_object()
  assign("scassist_compare_invalid_condition", obj, envir = .GlobalEnv)
  on.exit(rm("scassist_compare_invalid_condition", envir = .GlobalEnv), add = TRUE)

  expect_error(
    SCassist_compare_interactions(
      "scassist_compare_invalid_condition",
      group_by = "celltype",
      condition_by = "condition",
      condition_a = "missing",
      condition_b = "treated",
      run_llm = FALSE
    ),
    "condition_a"
  )
})

test_that("SCassist_compare_interactions errors when no shared cell groups exist", {
  obj <- scassist_compare_test_object(
    groups = c("A_only", "A_only", "A_only", "A_only", "B_only", "B_only", "B_only", "B_only")
  )
  assign("scassist_compare_no_shared_groups", obj, envir = .GlobalEnv)
  on.exit(rm("scassist_compare_no_shared_groups", envir = .GlobalEnv), add = TRUE)

  expect_error(
    SCassist_compare_interactions(
      "scassist_compare_no_shared_groups",
      group_by = "celltype",
      condition_by = "condition",
      condition_a = "control",
      condition_b = "treated",
      run_llm = FALSE
    ),
    "no shared cell groups"
  )
})

test_that("run_llm FALSE does not require an API key and preserves llm_response", {
  comparison_results <- list(
    llm_prompt = "Interpret this compact differential CellChat context.",
    llm_response = "previous response"
  )

  result <- SCassist:::SCassist_compare_interactions_add_llm_response(
    comparison_results = comparison_results,
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

  expect_true("llm_response" %in% names(result))
  expect_null(result$llm_response)
})

test_that("add_expression_support FALSE preserves expression support names with NULL values", {
  interactions <- data.frame(
    source = "Sender",
    target = "Receiver",
    interaction_name = "L1_R1",
    stringsAsFactors = FALSE
  )

  result <- SCassist:::SCassist_compare_interactions_expression_support_disabled(interactions)

  expect_true("expression_support" %in% names(result))
  expect_true("expression_support_subunits" %in% names(result))
  expect_true("differential_interactions_with_expression_support" %in% names(result))
  expect_null(result$expression_support)
  expect_null(result$expression_support_subunits)
})

test_that("expression-support helper works on a small synthetic expression matrix", {
  expression_data <- matrix(
    c(
      1, 1, 4, 4,
      2, 2, 2, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("L1", "R1"), paste0("cell", 1:4))
  )
  meta <- data.frame(
    celltype = rep("Sender", 4),
    condition = c("control", "control", "treated", "treated"),
    row.names = paste0("cell", 1:4),
    stringsAsFactors = FALSE
  )

  metrics <- SCassist:::SCassist_compare_interactions_single_gene_metrics(
    gene = "L1",
    expression_data = expression_data,
    meta = meta,
    group_by = "celltype",
    condition_by = "condition",
    cell_group = "Sender",
    condition_a = "control",
    condition_b = "treated",
    pseudocount = 1e-6
  )

  expect_gt(metrics$avg_expr_condition_b, metrics$avg_expr_condition_a)
  expect_gt(metrics$logfc_condition_b_vs_a, 0.25)
  expect_equal(metrics$pct_expr_condition_a, 1)
  expect_equal(metrics$pct_expr_condition_b, 1)
})

test_that("expression support labels strong supported interactions", {
  skip_if_not_installed("Seurat")

  counts <- matrix(
    c(
      1, 1, 0, 0, 8, 8, 0, 0,
      0, 0, 1, 1, 0, 0, 8, 8,
      1, 1, 1, 1, 1, 1, 1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("L1", "R1", "House"), paste0("cell", 1:8))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  obj$celltype <- rep(c("Sender", "Receiver", "Sender", "Receiver"), each = 2)
  obj$condition <- rep(c("control", "control", "treated", "treated"), each = 2)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)

  interactions <- data.frame(
    source = "Sender",
    target = "Receiver",
    ligand = "L1",
    receptor = "R1",
    interaction_name = "L1_R1",
    pathway_name = "TEST",
    score_condition_a = 0.2,
    score_condition_b = 1.2,
    delta_score = 1,
    pval_condition_a = 0.01,
    pval_condition_b = 0.01,
    status = "increased_in_condition_b",
    annotation = "Secreted Signaling",
    autocrine = FALSE,
    stringsAsFactors = FALSE
  )

  result <- SCassist:::SCassist_compare_interactions_compute_expression_support(
    seurat_object = obj,
    differential_interactions = interactions,
    group_by = "celltype",
    condition_by = "condition",
    condition_a = "control",
    condition_b = "treated",
    shared_cell_groups = c("Sender", "Receiver"),
    expression_support_assay = "RNA",
    expression_support_slot = "counts",
    expression_support_logfc_threshold = 0.25,
    expression_support_min_pct = 0.1,
    expression_support_pseudocount = 1e-6
  )

  expect_equal(result$expression_support$expression_support_label, "strong_expression_support")
  expect_true("expression_support_label" %in% names(result$differential_interactions_with_expression_support))
})

test_that("missing expression genes produce insufficient support without error", {
  skip_if_not_installed("Seurat")

  obj <- scassist_compare_test_object()
  interactions <- data.frame(
    source = "Sender",
    target = "Receiver",
    ligand = "MissingLigand",
    receptor = "MissingReceptor",
    interaction_name = "MissingLigand_MissingReceptor",
    pathway_name = "TEST",
    score_condition_a = 0,
    score_condition_b = 1,
    delta_score = 1,
    pval_condition_a = NA_real_,
    pval_condition_b = 0.01,
    status = "gained_in_condition_b",
    annotation = "Secreted Signaling",
    autocrine = FALSE,
    stringsAsFactors = FALSE
  )

  result <- SCassist:::SCassist_compare_interactions_compute_expression_support(
    seurat_object = obj,
    differential_interactions = interactions,
    group_by = "celltype",
    condition_by = "condition",
    condition_a = "control",
    condition_b = "treated",
    shared_cell_groups = c("Sender", "Receiver"),
    expression_support_assay = "RNA",
    expression_support_slot = "counts",
    expression_support_logfc_threshold = 0.25,
    expression_support_min_pct = 0.1,
    expression_support_pseudocount = 1e-6
  )

  expect_equal(result$expression_support$expression_support_label, "insufficient_expression_data")
  expect_match(result$expression_support$expression_support_note, "Missing ligand genes")
})

test_that("multi-subunit receptor with missing subunit is not strong support", {
  skip_if_not_installed("Seurat")

  counts <- matrix(
    c(
      1, 1, 0, 0, 8, 8, 0, 0,
      0, 0, 1, 1, 0, 0, 8, 8,
      1, 1, 1, 1, 1, 1, 1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("L1", "R1", "House"), paste0("cell", 1:8))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  obj$celltype <- rep(c("Sender", "Receiver", "Sender", "Receiver"), each = 2)
  obj$condition <- rep(c("control", "control", "treated", "treated"), each = 2)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)

  interactions <- data.frame(
    source = "Sender",
    target = "Receiver",
    ligand = "L1",
    receptor = "R1_R2",
    interaction_name = "L1_R1_R2",
    pathway_name = "TEST",
    score_condition_a = 0,
    score_condition_b = 1,
    delta_score = 1,
    pval_condition_a = NA_real_,
    pval_condition_b = 0.01,
    status = "gained_in_condition_b",
    annotation = "Secreted Signaling",
    autocrine = FALSE,
    stringsAsFactors = FALSE
  )

  result <- SCassist:::SCassist_compare_interactions_compute_expression_support(
    seurat_object = obj,
    differential_interactions = interactions,
    group_by = "celltype",
    condition_by = "condition",
    condition_a = "control",
    condition_b = "treated",
    shared_cell_groups = c("Sender", "Receiver"),
    expression_support_assay = "RNA",
    expression_support_slot = "counts",
    expression_support_logfc_threshold = 0.25,
    expression_support_min_pct = 0.1,
    expression_support_pseudocount = 1e-6
  )

  expect_false(result$expression_support$expression_support_label == "strong_expression_support")
  expect_true("R2" %in% result$expression_support_subunits$gene)
  expect_equal(
    result$expression_support_subunits$support_label[result$expression_support_subunits$gene == "R2"],
    "missing_gene"
  )
})

test_that("differential pathway summaries work on small synthetic tables", {
  pathways_a <- data.frame(
    pathway_name = c("CXCL", "MIF", "TGFb"),
    number_of_significant_interactions = c(2, 2, 1),
    total_communication_score = c(1, 2, 1.5),
    dominant_senders = c("Sender", "Receiver", "Sender"),
    dominant_receivers = c("Receiver", "Sender", "Receiver"),
    stringsAsFactors = FALSE
  )
  pathways_b <- data.frame(
    pathway_name = c("CXCL", "MIF", "IFN"),
    number_of_significant_interactions = c(3, 1, 2),
    total_communication_score = c(2, 0.5, 3),
    dominant_senders = c("Sender", "Receiver", "Sender"),
    dominant_receivers = c("Receiver", "Sender", "Receiver"),
    stringsAsFactors = FALSE
  )

  result <- SCassist:::SCassist_compare_interactions_compare_pathways(pathways_a, pathways_b)
  status <- stats::setNames(result$status, result$pathway_name)

  expect_equal(status[["CXCL"]], "increased_in_condition_b")
  expect_equal(status[["MIF"]], "decreased_in_condition_b")
  expect_equal(status[["TGFb"]], "lost_in_condition_b")
  expect_equal(status[["IFN"]], "gained_in_condition_b")
})

test_that("differential interaction summaries work on small synthetic tables", {
  interactions_a <- data.frame(
    source = c("Sender", "Sender", "Receiver"),
    target = c("Receiver", "Receiver", "Receiver"),
    ligand = c("L1", "L2", "L3"),
    receptor = c("R1", "R2", "R3"),
    interaction_name = c("L1_R1", "L2_R2", "L3_R3"),
    pathway_name = c("CXCL", "MIF", "TGFb"),
    communication_score = c(1, 2, 0.5),
    pval = c(0.01, 0.02, 0.03),
    annotation = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
    autocrine = c(FALSE, FALSE, TRUE),
    rank_global = 1:3,
    stringsAsFactors = FALSE
  )
  interactions_b <- data.frame(
    source = c("Sender", "Sender", "Receiver"),
    target = c("Receiver", "Receiver", "Receiver"),
    ligand = c("L1", "L4", "L3"),
    receptor = c("R1", "R4", "R3"),
    interaction_name = c("L1_R1", "L4_R4", "L3_R3"),
    pathway_name = c("CXCL", "IFN", "TGFb"),
    communication_score = c(2, 3, 0.25),
    pval = c(0.01, 0.02, 0.03),
    annotation = c("Secreted Signaling", "Secreted Signaling", "Cell-Cell Contact"),
    autocrine = c(FALSE, FALSE, TRUE),
    rank_global = 1:3,
    stringsAsFactors = FALSE
  )

  result <- SCassist:::SCassist_compare_interactions_compare_ligand_receptors(interactions_a, interactions_b)
  status <- stats::setNames(result$status, result$interaction_name)

  expect_equal(status[["L1_R1"]], "increased_in_condition_b")
  expect_equal(status[["L2_R2"]], "lost_in_condition_b")
  expect_equal(status[["L4_R4"]], "gained_in_condition_b")
  expect_equal(status[["L3_R3"]], "decreased_in_condition_b")
  expect_true(result$autocrine[result$interaction_name == "L3_R3"])
})

test_that("differential LLM context and prompt contain comparison caveats", {
  metadata <- list(
    condition_a = "control",
    condition_b = "treated",
    shared_cell_groups = c("Sender", "Receiver"),
    condition_a_only_cell_groups = character(),
    condition_b_only_cell_groups = "Responder"
  )
  global_comparison <- list(delta_interaction_count = 3)
  differential_pathways <- SCassist:::SCassist_compare_interactions_empty_pathway_comparison()
  differential_interactions <- SCassist:::SCassist_compare_interactions_empty_lr_comparison()
  differential_cell_roles <- SCassist:::SCassist_compare_interactions_empty_role_comparison()
  gained_lost_summary <- SCassist:::SCassist_compare_interactions_empty_gained_lost_summary()
  expression_support <- data.frame(
    source = "Sender",
    target = "Receiver",
    interaction_name = "L1_R1",
    pathway_name = "TEST",
    status = "increased_in_condition_b",
    score_condition_a = 0.2,
    score_condition_b = 1.2,
    delta_score = 1,
    ligand = "L1",
    receptor = "R1",
    ligand_genes = "L1",
    receptor_genes = "R1",
    missing_ligand_genes = NA_character_,
    missing_receptor_genes = NA_character_,
    ligand_avg_expr_condition_a = 1,
    ligand_avg_expr_condition_b = 4,
    ligand_logfc_condition_b_vs_a = 2,
    ligand_pct_expr_condition_a = 1,
    ligand_pct_expr_condition_b = 1,
    ligand_delta_pct_expr = 0,
    receptor_avg_expr_condition_a = 1,
    receptor_avg_expr_condition_b = 4,
    receptor_logfc_condition_b_vs_a = 2,
    receptor_pct_expr_condition_a = 1,
    receptor_pct_expr_condition_b = 1,
    receptor_delta_pct_expr = 0,
    expected_direction = "increase_in_condition_b",
    ligand_expression_support = "supports_expected_direction",
    receptor_expression_support = "supports_expected_direction",
    expression_support_label = "strong_expression_support",
    expression_support_note = "synthetic strong support",
    stringsAsFactors = FALSE
  )

  llm_context <- SCassist:::SCassist_compare_interactions_build_llm_context(
    metadata = metadata,
    global_comparison = global_comparison,
    differential_pathways = differential_pathways,
    differential_interactions = differential_interactions,
    differential_cell_roles = differential_cell_roles,
    gained_lost_summary = gained_lost_summary,
    expression_support = expression_support,
    expression_support_subunits = SCassist:::SCassist_compare_interactions_empty_expression_support_subunits(),
    top_n_interactions = 10,
    top_n_pathways = 10,
    experimental_context = "synthetic test"
  )
  llm_prompt <- SCassist:::SCassist_compare_interactions_build_llm_prompt(llm_context)

  expect_true("interpretation_caveats" %in% names(llm_context))
  expect_true("expression_support_summary" %in% names(llm_context))
  expect_match(llm_prompt, "pairwise differential")
  expect_match(llm_prompt, "shared cell groups only")
  expect_match(llm_prompt, "not formal DEG-supported")
  expect_match(llm_prompt, "Expression support is based on transcript-level")
})

test_that("CellChat-dependent comparison setup smoke test skips cleanly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("CellChat")

  obj <- scassist_compare_test_object()
  condition_info <- SCassist:::SCassist_compare_interactions_get_condition_groups(
    seurat_object = obj,
    group_by = "celltype",
    condition_by = "condition",
    condition_a = "control",
    condition_b = "treated"
  )
  db_info <- SCassist:::SCassist_interactions_get_cellchat_db("human", "default")

  expect_equal(condition_info$shared_cell_groups, c("Receiver", "Sender"))
  expect_true("interaction" %in% names(db_info$database))
})
