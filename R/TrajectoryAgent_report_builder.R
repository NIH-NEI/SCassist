# JSON, LLM-context, and markdown report builder for TrajectoryAgent.

ta_cautious_interpretation_rules <- function() {
  c(
    "Use only the structured evidence package.",
    "Do not infer causality.",
    "Do not claim Monocle3 proves differentiation.",
    "Do not claim lineage commitment unless branch evidence is strong.",
    "Use cautious phrases such as consistent with, suggesting, may reflect, and supports but does not prove."
  )
}

ta_build_llm_context <- function(evidence) {
  list(
    trajectory_structure = evidence$trajectory_structure,
    root_selection_evidence = evidence$root_selection_evidence,
    pseudotime_qc = evidence$pseudotime_qc,
    branch_summary = evidence$branch_summary,
    dynamic_gene_summary = evidence$dynamic_gene_summary,
    gene_module_summary = evidence$gene_module_summary,
    metadata_associations = evidence$metadata_associations,
    confidence_flags = evidence$confidence_flags,
    confidence_score = evidence$confidence_score,
    confidence_label = evidence$confidence_label,
    cautious_interpretation_rules = ta_cautious_interpretation_rules()
  )
}

ta_deterministic_interpretation <- function(evidence) {
  executive <- character()
  if (length(evidence$root_selection_evidence$root_pr_nodes) > 0 &&
    evidence$pseudotime_qc$finite_pseudotime_cells > 0) {
    executive <- c(executive, "Cells are ordered along a trajectory consistent with the supplied root evidence.")
  } else {
    executive <- c(executive, "Trajectory direction is low confidence because no defensible root was selected.")
  }
  executive <- c(
    executive,
    paste0(
      "The principal graph contains ",
      evidence$trajectory_structure$num_branch_points,
      " branch point(s), ",
      evidence$trajectory_structure$num_leaf_nodes,
      " leaf node(s), and ",
      evidence$trajectory_structure$num_partitions,
      " partition(s)."
    )
  )
  if (evidence$pseudotime_qc$infinite_pseudotime_cells > 0) {
    executive <- c(
      executive,
      paste0(evidence$pseudotime_qc$infinite_pseudotime_cells, " cells have infinite pseudotime and should be treated as a QC signal.")
    )
  }

  branches <- if (length(evidence$branch_summary) == 0) {
    "No interpretable branch point was detected in the principal graph."
  } else {
    vapply(evidence$branch_summary, function(branch) {
      paste0("Branch node ", branch$branch_node_id, " has confidence label ", branch$confidence, "; branch-associated states should be interpreted cautiously.")
    }, character(1))
  }
  dynamic <- if (evidence$dynamic_gene_summary$num_significant_genes > 0) {
    paste0(
      evidence$dynamic_gene_summary$num_significant_genes,
      " genes vary over the trajectory by graph_test, suggesting dynamic expression programs along the graph."
    )
  } else {
    "No trajectory-variable genes passed the configured graph_test thresholds."
  }
  caveats <- c(
    "Monocle3 ordering supports, but does not prove, developmental progression.",
    "Pseudotime direction depends on the selected root evidence.",
    "Branch-associated modules may reflect state differences, technical effects, or sampling imbalance."
  )
  if (length(evidence$metadata_associations$possible_batch_or_donor_confounding) > 0) {
    caveats <- c(caveats, "Batch or donor metadata is associated with pseudotime and should be treated as a potential confounder.")
  }
  followup <- c(
    "Validate root choice using independent time, marker, or experimental prior information.",
    "Inspect UMAP, principal graph, and pseudotime overlays before making biological claims.",
    "Test key dynamic genes or modules in independent data or orthogonal assays."
  )
  list(
    executive_summary_candidate = executive,
    branch_interpretation_candidates = branches,
    dynamic_program_candidates = dynamic,
    caveats = caveats,
    followup_suggestions = followup
  )
}

ta_write_trajectory_reports <- function(evidence, output_dir) {
  ta_validate_evidence_schema(evidence)
  llm_context <- ta_build_llm_context(evidence)
  ta_validate_llm_context(llm_context)
  ta_write_json(evidence, file.path(output_dir, "trajectory_evidence.json"))
  ta_write_json(llm_context, file.path(output_dir, "llm_context.json"))

  interp <- evidence$llm_interpretation_input
  lines <- c(
    "# SCAssist TrajectoryAgent Report",
    "",
    "## 1. Executive summary",
    paste0("- ", interp$executive_summary_candidate),
    "",
    "## 2. Monocle3 run configuration",
    "```json",
    jsonlite::toJSON(evidence$run_config, pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null"),
    "```",
    "",
    "## 3. Root selection evidence",
    paste0("- Method: ", evidence$root_selection_evidence$method),
    paste0("- Root nodes: ", paste(evidence$root_selection_evidence$root_pr_nodes, collapse = ", ")),
    paste0("- Partitions without roots: ", paste(evidence$root_selection_evidence$partitions_without_roots, collapse = ", ")),
    paste0("- Confidence: ", evidence$root_selection_evidence$confidence),
    "",
    "## 4. Pseudotime QC",
    paste0("- Finite pseudotime cells: ", evidence$pseudotime_qc$finite_pseudotime_cells),
    paste0("- Infinite pseudotime cells: ", evidence$pseudotime_qc$infinite_pseudotime_cells),
    paste0("- Fraction infinite pseudotime: ", round(evidence$pseudotime_qc$fraction_infinite_pseudotime, 3)),
    "",
    "## 5. Trajectory structure",
    paste0("- Cells: ", evidence$trajectory_structure$num_cells),
    paste0("- Genes: ", evidence$trajectory_structure$num_genes),
    paste0("- Clusters: ", evidence$trajectory_structure$num_clusters),
    paste0("- Partitions: ", evidence$trajectory_structure$num_partitions),
    "",
    "## 6. Partition summary",
    paste0("- Cells per partition: ", paste(names(evidence$trajectory_structure$cells_per_partition), evidence$trajectory_structure$cells_per_partition, sep = "=", collapse = ", ")),
    "",
    "## 7. Branch and leaf summary",
    paste0("- Branch points: ", paste(evidence$trajectory_structure$branch_point_ids, collapse = ", ")),
    paste0("- Leaf nodes: ", paste(evidence$trajectory_structure$leaf_node_ids, collapse = ", ")),
    paste0("- ", interp$branch_interpretation_candidates),
    "",
    "## 8. Dynamic trajectory genes",
    paste0("- ", interp$dynamic_program_candidates),
    "",
    "## 9. Gene modules",
    paste0("- Detected modules: ", evidence$gene_module_summary$num_modules),
    "- Branch-associated modules may reflect coordinated trajectory-variable programs when supported by module activity summaries.",
    "",
    "## 10. Metadata associations and possible confounding",
    paste0("- Possible batch/donor confounding: ", paste(evidence$metadata_associations$possible_batch_or_donor_confounding, collapse = ", ")),
    "",
    "## 11. Confidence assessment",
    paste0("- Score: ", evidence$confidence_score),
    paste0("- Label: ", evidence$confidence_label),
    paste0("- Flags: ", paste(evidence$confidence_flags, collapse = ", ")),
    "",
    "## 12. Cautionary notes",
    paste0("- ", interp$caveats),
    "",
    "## 13. Suggested follow-up analyses",
    paste0("- ", interp$followup_suggestions),
    ""
  )
  writeLines(lines, file.path(output_dir, "trajectory_report.md"), useBytes = TRUE)
  invisible(llm_context)
}
