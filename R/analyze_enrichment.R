#' @title Analyzes KEGG and GO Pathway Enrichment Results and Provides Insights from a Language Model
#'
#' @description This function analyzes KEGG and/or GO pathway enrichment results obtained from a list of genes
#' provided by the user, and then uses a language model 
#' to generate insights based on the results. 
#' 
#' @author Vijay Nagarajan, PhD, NEI/NIH
#'
#' @details This function was written with assistance from Google's Gemini and Meta's Llama3.
#'
#' @param organism The organism of interest. Supported organisms are:
#'   "human", "mouse", "rat", "zebrafish", "fly", "yeast", "worm", "arabidopsis",
#'   "chicken". Default is "human".
#' @param model_G Character string specifying the Gemini model to use for
#'  analysis. Default is "gemini-1.5-flash-latest".
#' @param model_O Character string specifying the Ollama model to use for
#'  analysis. Default is "llama3".
#' @param model_C Character string specifying the OpenAI model to use for
#'  analysis. Default is "gpt-4o-mini".
#' @param markers A data frame containing the marker genes for a cluster,
#'   with columns for 'p_val', and 'avg_log2FC'. **Gene names should be present
#'   as row names of the data frame.** This data frame can be the output of the `FindMarkers`
#'   function from Seurat.
#' @param pvalue The p-value threshold for selecting significant genes. Default is 0.05.
#' @param log2FC The absolute value of the average log2 fold change 
#'   threshold for selecting significant genes. Default is 1.
#' @param experimental_design A brief description of the experimental design, e.g.,
#'   "Single-cell RNA sequencing of immune cells from patients with cancer". Default is "Single-cell RNA sequencing".
#' @param llm_server The LLM server to use. Options are "google" or "ollama" or "openai". Default is "google".
#' @param output_file The name of the text file to save the network data to.
#'   Default is "L_SCassist_summary_network_data.txt".
#' @param html_file The name of the HTML file to save the network visualization to.
#'   Default is "L_SCassist_summary_network_visualization.html".
#' @param do_kegg Whether to perform KEGG enrichment analysis. Default is TRUE.
#' @param do_go Whether to perform GO enrichment analysis. Default is TRUE.
#' @param go_ont The ontology to use for GO enrichment. Options are "BP", "CC", "MF", or "ALL". Default is "BP".
#' @param temperature The temperature parameter for the language model. Higher
#' values lead to more creative and unpredictable outputs. Default is 0.
#' @param max_output_tokens The maximum number of tokens the language model can generate.
#' Default is 10024.
#' @param api_key_file The path to the file containing the API key.
#' @param model_params A list of parameters to be passed to the `rollama::query` function.
#'   This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.
#' 
#' @usage
#' # Using Google's Gemini Pro:
#' SCassist_analyze_enrichment(llm_server="google", organism = "human", markers, 
#'          pvalue = 0.05, log2FC = 1, 
#'          experimental_design = "Single-cell RNA sequencing",
#'          temperature = 0,
#'          max_output_tokens = 10048,
#'          model_G = "gemini-1.5-flash-latest",
#'          model_O = "llama3",
#'          model_C = "gpt-4o-mini",
#'          api_key_file = "api_keys.txt",
#'          output_file = "G_SCassist_summary_network_data.txt", 
#'          html_file = "G_SCassist_summary_network_visualization.html",
#'          do_kegg = TRUE, do_go = TRUE, go_ont = "BP",
#'          model_params = list(seed = 42, temperature = 0, num_gpu = 0))
#' 
#' @return The function prints a response from a language model that provides insights
#'   on the enriched pathways, their relationships, the biological system, and potential
#'   key genes or targets. It also saves the enrichment results to a text file
#'   named "enrichment_results.txt" in the current working directory.
#'
#' @keywords KEGG, GO, pathway enrichment, language model, single-cell RNA sequencing
#'
#'
#' @import BiocManager
#' @import org.Hs.eg.db
#' @import rollama
#' @import visNetwork
#' @import jsonlite
#' @importFrom dplyr %>%
#' @importFrom utils capture.output head
#' @importFrom data.table as.data.table
#' @importFrom clusterProfiler enrichKEGG bitr enrichGO
#' @importFrom httr POST content content_type_json
#' @importFrom jsonlite fromJSON
#' @export
#' 

SCassist_analyze_enrichment <- function(llm_server="google", organism = "human", markers = "", 
                                        pvalue = 0.05, log2FC = 1, 
                                        experimental_design = "Single-cell RNA sequencing",
                                        temperature = 0,
                                        max_output_tokens = 10048,
                                        model_G = "gemini-1.5-flash-latest",
                                        model_O = "llama3",
                                        model_C = "gpt-4o-mini",
                                        api_key_file = "api_keys.txt",
                                        output_file = "G_SCassist_summary_network_data.txt", 
                                        html_file = "G_SCassist_summary_network_visualization.html",
                                        do_kegg = TRUE, do_go = TRUE, go_ont = "BP",
                                        model_params = list(seed = 42, temperature = 0, num_gpu = 0)
) {
  
  if (llm_server == "google") {
    return(SCassist_analyze_enrichment_G(organism = organism, markers = markers, 
                                         pvalue = pvalue, log2FC = log2FC, 
                                         experimental_design = experimental_design,
                                         temperature = temperature,
                                         max_output_tokens = max_output_tokens,
                                         model = model_G,
                                         api_key_file = api_key_file,
                                         output_file = output_file, 
                                         html_file = html_file,
                                         do_kegg = do_kegg, do_go = do_go, go_ont = go_ont))
  } else if (llm_server == "ollama") {
    return(SCassist_analyze_enrichment_L(organism = organism, markers = markers, 
                                         pvalue = pvalue, log2FC = log2FC, 
                                         experimental_design = experimental_design,
                                         output_file = output_file, 
                                         html_file = html_file,
                                         model = model_O,
                                         do_kegg = do_kegg, do_go = do_go, go_ont = go_ont,
                                         model_params = model_params
    ))
  } else if (llm_server == "openai") {
    return(SCassist_analyze_enrichment_C(organism = organism, markers = markers, 
                                         pvalue = pvalue, log2FC = log2FC, 
                                         experimental_design = experimental_design,
                                         temperature = temperature,
                                         max_output_tokens = max_output_tokens,
                                         model = model_C,
                                         api_key_file = api_key_file,
                                         output_file = output_file, 
                                         html_file = html_file,
                                         model_params = model_params,
                                         do_kegg = do_kegg, do_go = do_go, go_ont = go_ont
    ))
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama'.")
  }
}

SCassist_analyze_enrichment_G <- function(organism = "human", markers="", 
                                          pvalue = 0.05, log2FC = 1, 
                                          experimental_design = "Single-cell RNA sequencing",
                                          temperature = 0,
                                          max_output_tokens = 10048,
                                          model = "gemini-1.5-flash-latest",
                                          api_key_file = "api_keys.txt",
                                          output_file = "G_SCassist_summary_network_data.txt", 
                                          html_file = "G_SCassist_summary_network_visualization.html",
                                          do_kegg = TRUE, do_go = TRUE, go_ont = "BP") {
  
  model_query <- paste0(model, ":generateContent")
  
  Description <- Gene <- SYMBOL <- geneID <- V <- subcategory <- . <- NULL
  
  # Read the API key from the specified file
  api_key <- readLines(api_key_file)
  
  # Suppress the warning message
  suppressWarnings({
    # Function to check and install the database if needed
    install_organism_db <- function(organism) {
      organism_db <- switch(organism,
                            "human" = "org.Hs.eg.db",
                            "mouse" = "org.Mm.eg.db",
                            "rat" = "org.Rn.eg.db",
                            "zebrafish" = "org.Dr.eg.db",
                            "fly" = "org.Dm.eg.db",
                            "yeast" = "org.Sc.sgd.db",
                            "worm" = "org.Ce.eg.db",
                            "arabidopsis" = "org.At.tair.db",
                            "chicken" = "org.Gg.eg.db",
                            stop("Unsupported organism. Please choose from the list.")
      )
      
      if (!requireNamespace(organism_db, quietly = TRUE)) {
        cat(paste("The", organism_db, "database is not installed. Installing...\n"))
        BiocManager::install(organism_db)
      }
      
      # Use :: to access the library function from the specific package
      library(organism_db, character.only = TRUE)
      return(organism_db)
    }
    
    # Convert organism name to a standard code used by clusterProfiler
    organism_code <- switch(organism,
                            "human" = "hsa",
                            "mouse" = "mmu",
                            "rat" = "rno",
                            "zebrafish" = "dre",
                            "fly" = "dme",
                            "yeast" = "sce",
                            "worm" = "cel",
                            "arabidopsis" = "ath",
                            "chicken" = "gga",
                            stop("Unsupported organism. Please choose from the list.")
    )
    
    # 1. Check and install the database if needed
    organism_db <- install_organism_db(organism)
    
    # 2. Use :: to access the functions from the specific packages
    markers <- markers %>%
      dplyr::filter(markers[["p_val"]] < pvalue & (abs(markers[["avg_log2FC"]]) > log2FC))
    
    # 3. Extract gene names as a list
    gene_list <- rownames(markers) 
    gene_list <- sub("\\..*", "", gene_list)
    
    # 4. Determine input type (SYMBOL or ENTREZID)
    input_type <- ifelse(is.character(gene_list[[1]]), "SYMBOL", "ENTREZID")
    
    # 5. Create a table with gene names and IDs based on input type
    
    if (input_type == "ENTREZID") {
      gene_mapping <- suppressMessages(clusterProfiler::bitr(gene_list, fromType = "ENTREZID", toType = c("SYMBOL", "ENTREZID"), OrgDb = organism_db))
    } else {
      gene_mapping <- suppressMessages(clusterProfiler::bitr(gene_list, fromType = "SYMBOL", toType = c("SYMBOL", "ENTREZID"), OrgDb = organism_db))
    }
    
    # KEGG Enrichment
    if (do_kegg) {
      
      writeLines("\nThe KEGG analysis is presented below.
-------------------------------------------\n")
      
      # 6. Perform pathway enrichment analysis using ENTREZID
      kegg_enrichment <- clusterProfiler::enrichKEGG(gene = gene_mapping$ENTREZID, organism = organism_code)
      
      # 7. Get the output of kegg enrichment analysis into a separate table
      kegg_enrichment_results <- kegg_enrichment@result

      # 13. Save enrichment results to a text file
      write.table(kegg_enrichment_results, file = "G_kegg_enrichment_results_original.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The ClusterProfiler KEGG pathways result is saved as a text file in the current working directory. \n")
      
      # 8. Filter for pvalue < 0.05
      kegg_enrichment_results <- kegg_enrichment_results %>%
        dplyr::filter(pvalue < 0.05)
      
      # 9. Map gene names to gene IDs
      kegg_enrichment_results <- kegg_enrichment_results %>%
        dplyr::mutate(geneID = strsplit(geneID, "/")) %>% 
        tidyr::unnest(geneID) %>%
        # Convert geneID to character for joining
        dplyr::mutate(geneID = as.character(geneID)) %>%
        # Join with gene_mapping based on geneID and ENTREZID
        dplyr::left_join(gene_mapping, by = c("geneID" = "ENTREZID")) %>%
        # Replace geneID with SYMBOL if available
        dplyr::mutate(Gene = ifelse(!is.na(SYMBOL), SYMBOL, geneID)) %>% 
        # Group by subcategory, Description, and pvalue
        dplyr::group_by(subcategory, Description, pvalue) %>%
        # Summarize the results, keeping only unique combinations
        dplyr::summarize(., geneID = paste(unique(geneID), collapse = ", "),
                         Gene = paste(unique(Gene), collapse = ", "), .groups = "drop") 
      
      # Use the square bracket notation
      input_kegg_enrichment <- kegg_enrichment_results[, c("Description", "pvalue", "Gene")]
      
      # 13. Save enrichment results to a text file
      write.table(kegg_enrichment_results, file = "G_kegg_enrichment_results_genename.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The list of significantly enriched KEGG pathways is also saved as a text file in the current working directory. \n")

      json_kegg_enrichment <- jsonlite::toJSON(input_kegg_enrichment)
      
      input_text_kegg <- paste0("
This below Data is a list of KEGG pathway enrichment results for a set of differentially expressed genes obtained from a single cell experiment involving ",experimental_design,".\n 
Data:
```json \n",
                                json_kegg_enrichment,
                                "```\n
Analyze all of the pathways from my above Data and provide insights in a structured format. Include:
                           
        1. **Significant Pathways:** Analyze all the pathways in my list, in the context of the system involving ",experimental_design," and summarize common themes, in 2 paragraph, with a total of 10 lines.
        
        2. **Regulators:** Include potential involvment of any transcription factors from the given list, in the context of the system involving ",experimental_design," in a 5 line paragraph.
        
        3. **Key Genes or Potential Targets:** Suggest which genes from the enriched pathways, in my data input, might be important to the system or a target based on their potential impact on the system involving", experimental_design," in a 10 line paragraph.
        
        Do not provide any Further Investigation or comments, not asked for.")
      
      # 10. Send the request to the Gemini API
      response_kegg <- POST(
        url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
        query = list(key = api_key), # Access api_key from api_keys.R
        content_type_json(),
        encode = "json",
        body = list(
          contents = list(
            parts = list(
              list(text = input_text_kegg)
            )),
          generationConfig = list(
            temperature = temperature,
            maxOutputTokens = max_output_tokens,
            seed = 123456  # For reproducible runs
          )
        )
      )
      
      if(response_kegg$status_code>200) {
        stop(paste("Error - ", content(response_kegg)$error$message))
      }
      
      # 11. Extract the Gemini model's response
      candidates_kegg <- httr::content(response_kegg)$candidates
      outputs_kegg <- unlist(lapply(candidates_kegg, function(candidate) candidate$content$parts))
      
      # 12. Print the LLM's response
      cat("KEGG Enrichment Insights:\n")
      cat(outputs_kegg)
      
      # 14. Return the raw output from the Gemini Pro language model
      input_for_network_kegg <- paste0("Extract important named entities and relationships between key genes and 
                                potential concepts, not cell types, and format it in three columns as Gene, 
                                Interaction, Concept, one gene - one interaction per row. 
                                the gene column should only contain one gene. Include only concepts related to biological pathways.
                                the concept column should only contain ONE gene or ONE concept. 
                                Add potential gene-gene interactions as well, for the core system. 
                                Do not provide any concept - concept associations.
                                EXAMPLE OUTPUT:
                                | Gene1 | Involved in | Metabolism |
                                | Gene3 | Interacts with | Gene5 |\n",
                                       outputs_kegg)
      
      # 10. Send the request to the Gemini API
      response2_kegg <- POST(
        url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
        query = list(key = api_key), # Access api_key from api_keys.R
        content_type_json(),
        encode = "json",
        body = list(
          contents = list(
            parts = list(
              list(text = input_for_network_kegg)
            )),
          generationConfig = list(
            temperature = temperature,
            maxOutputTokens = max_output_tokens,
            seed = 123456  # For reproducible runs
          )
        )
      )
      
      if(response2_kegg$status_code>200) {
        stop(paste("Error - ", content(response2_kegg)$error$message))
      }
      
      # 11. Extract the Gemini model's response
      candidates2_kegg <- httr::content(response2_kegg)$candidates
      outputs2_kegg <- unlist(lapply(candidates2_kegg, function(candidates2_kegg) candidates2_kegg$content$parts))
      
      # 12. Print the LLM's response
      
      # Extract and format the network data from the LLM response
      network_section_kegg <- strsplit(outputs2_kegg, "\n")[[1]]
      network_data_kegg <- network_section_kegg[5:length(network_section_kegg)]
      network_data_kegg <- strsplit(network_data_kegg, "\\|\\s*")

      #cat(paste(network_data_kegg, collapse = "\n"))
      #str(network_data_kegg) 
      #class(network_data_kegg)
      
      # Create a data frame from the network data
      network_df_kegg <- data.frame(
        Gene = sapply(network_data_kegg, "[", 2),
        Interaction = sapply(network_data_kegg, "[", 3),
        Concept = sapply(network_data_kegg, "[", 4)
      )
      
      # Clean up the data frame
      network_df_kegg <- data.frame(
        lapply(network_df_kegg, function(x) gsub("\\s+", " ", trimws(x)))
      )
    }
    
    # GO Enrichment
    if (do_go) {
      
      writeLines("\nThe GO analysis is provided below.
----------------------------------------\n")
      
      # 6. Perform pathway enrichment analysis using ENTREZID
      go_enrichment <- clusterProfiler::enrichGO(gene = gene_mapping$ENTREZID, OrgDb = organism_db, ont=go_ont)
      
      # 7. Get the output of kegg enrichment analysis into a separate table
      go_enrichment_results <- go_enrichment@result
      
      # 13. Save enrichment results to a text file
      write.table(go_enrichment_results, file = "G_go_enrichment_results_original.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The ClusterProfiler GO enrichment result is saved as a text file in the current working directory. \n")
      
      # 8. Filter for pvalue < 0.05
      go_enrichment_results <- go_enrichment_results %>%
        dplyr::filter(pvalue < 0.05)
      
      # 9. Map gene names to gene IDs
      go_enrichment_results <- go_enrichment_results %>%
        dplyr::mutate(geneID = strsplit(geneID, "/")) %>% 
        tidyr::unnest(geneID) %>%
        # Convert geneID to character for joining
        dplyr::mutate(geneID = as.character(geneID)) %>%
        # Join with gene_mapping based on geneID and ENTREZID
        dplyr::left_join(gene_mapping, by = c("geneID" = "ENTREZID")) %>%
        # Replace geneID with SYMBOL if available
        dplyr::mutate(Gene = ifelse(!is.na(SYMBOL), SYMBOL, geneID)) %>% 
        # Group by subcategory, Description, and pvalue
        dplyr::group_by(Description, pvalue) %>%
        # Summarize the results, keeping only unique combinations
        dplyr::summarize(., geneID = paste(unique(geneID), collapse = ", "),
                         Gene = paste(unique(Gene), collapse = ", "), .groups = "drop") 
      
      # Alternatively, use the square bracket notation
      input_go_enrichment <- go_enrichment_results[, c("Description", "pvalue", "Gene")]

      # 13. Save enrichment results to a text file
      write.table(go_enrichment_results, file = "G_go_enrichment_results_genename.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The list of significantly enriched GO terms is also saved as a text file in the current working directory. \n")
      
      json_go_enrichment <- jsonlite::toJSON(input_go_enrichment)
      
      input_text_go <- paste0("
This below Data is a list of GO enrichment results for a set of differentially expressed genes obtained from a single cell experiment involving ",experimental_design,".\n 
Data:
```json \n",
                              json_go_enrichment,
                              "```\n
Analyze all of the enriched concepts from my above Data and provide insights in a structured format. Include:
                           
        1. **Significant Concepts:** Analyze all the concepts in my list, in the context of the system involving ",experimental_design," and summarize common themes, in 2 paragraph, with a total of no more than 10 lines.
        
        2. **Regulators:** Include potential involvment of any transcription factors from the given list, in the context of the system involving ",experimental_design," in a 5 line paragraph. Restrict the number of genes to less than 10.
        
        3. **Key Genes or Potential Targets:** Suggest which genes from the enriched ontologies, in my data input, might be important to the system or a target based on their potential impact on the system involving", experimental_design," in a 10 line paragraph.
        
        Do not provide any Further Investigation or comments, not asked for.")
      
      # 10. Send the request to the Gemini API
      response_go <- POST(
        url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
        query = list(key = api_key),
        content_type_json(),
        encode = "json",
        body = list(
          contents = list(
            parts = list(
              list(text = input_text_go)
            )),
          generationConfig = list(
            temperature = temperature,
            maxOutputTokens = max_output_tokens,
            seed = 123456  # For reproducible runs
          )
        )
      )
      
      if(response_go$status_code>200) {
        stop(paste("Error - ", content(response_go)$error$message))
      }
      
      # 11. Extract the Gemini model's response
      candidates_go <- httr::content(response_go)$candidates
      outputs_go <- unlist(lapply(candidates_go, function(candidate) candidate$content$parts))
      
      # 12. Print the LLM's response
      cat("GO Enrichment Insights:\n")
      cat(outputs_go)
      
      # 14. Return the raw output from the Gemini Pro language model
      input_for_network_go <- paste0("Extract important named entities and relationships between key genes and 
                                potential concepts, not cell types, and format it in three columns as Gene, 
                                Interaction, Concept, one gene - one interaction per row. 
                                the gene column should only contain one gene. Include only concepts related to GO terms.
                                the concept column should only contain ONE gene or ONE concept. 
                                Add potential gene-gene interactions as well, for the core system. 
                                Do not provide any concept - concept associations.
                                EXAMPLE OUTPUT:
                                | Gene1 | Involved in | Metabolism |
                                | Gene3 | Interacts with | Gene5 |\n",
                                     outputs_go)
      
      # 10. Send the request to the Gemini API
      response2_go <- POST(
        url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
        query = list(key = api_key), # Access api_key from api_keys.R
        content_type_json(),
        encode = "json",
        body = list(
          contents = list(
            parts = list(
              list(text = input_for_network_go)
            )),
          generationConfig = list(
            temperature = temperature,
            maxOutputTokens = max_output_tokens,
            seed = 123456  # For reproducible runs
          )
        )
      )
      
      if(response2_go$status_code>200) {
        stop(paste("Error - ", content(response2_go)$error$message))
      }
      
      # 11. Extract the Gemini model's response
      candidates2_go <- httr::content(response2_go)$candidates
      outputs2_go <- unlist(lapply(candidates2_go, function(candidates2_go) candidates2_go$content$parts))
      
      # 12. Print the LLM's response
      
      # Extract and format the network data from the LLM response
      network_section_go <- strsplit(outputs2_go, "\n")[[1]]
      network_data_go <- network_section_go[5:length(network_section_go)]
      network_data_go <- strsplit(network_data_go, "\\|\\s*")
      
      # Create a data frame from the network data
      network_df_go <- data.frame(
        Gene = sapply(network_data_go, "[", 2),
        Interaction = sapply(network_data_go, "[", 3),
        Concept = sapply(network_data_go, "[", 4)
      )
      
      # Clean up the data frame
      network_df_go <- data.frame(
        lapply(network_df_go, function(x) gsub("\\s+", " ", trimws(x)))
      )
    }
    
    # Combine Networks
    if (do_kegg && do_go) {
      network_df_combined <- rbind(network_df_kegg, network_df_go)
      
      # Save the combined network data to a tab-delimited text file
      write.table(network_df_combined, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create an interactive network visualization
      network_df_combined <- data.frame(
        source = network_df_combined$Gene,
        target = network_df_combined$Concept,
        interaction = network_df_combined$Interaction
      )
      
      # Create the graph object
      network_graph_combined <- igraph::graph_from_data_frame(network_df_combined, directed = TRUE)
      
      # Create the visNetwork object
      vis_network_combined <- visNetwork(nodes = data.frame(id = igraph::V(network_graph_combined)$name, 
                                                            label = igraph::V(network_graph_combined)$name), 
                                         edges = data.frame(from = network_df_combined$source, 
                                                            to = network_df_combined$target, 
                                                            label = network_df_combined$interaction))
      
      # Customize the visualization
      vis_network_combined <- vis_network_combined %>%
        visEdges(arrows = "to") %>% 
        visNodes(size = 20,  # Adjust node size as needed
                 font = list(size = 24)) %>%  # Increase node label size
        visOptions(highlightNearest = TRUE)  # Highlight nodes when hovering
      
      # Save the visualization as an HTML file
      visSave(vis_network_combined, file = html_file)
      
      # Combine Gemini Pro responses for overall summary
      combined_response <- paste0(outputs_kegg, "\n\n", outputs_go)
      
      # Generate overall summary prompt
      input_text_overall <- paste0("
This is a combined summary of KEGG and GO enrichment results for a set of differentially expressed genes obtained from a single cell experiment involving ",experimental_design,".\n 
Please provide a comprehensive summary based on the following information:
", combined_response, "\n
Do not provide any Further Investigation or comments, not asked for.")
      
      # Send the request to the Gemini API
      response_overall <- POST(
        url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
        query = list(key = api_key),
        content_type_json(),
        encode = "json",
        body = list(
          contents = list(
            parts = list(
              list(text = input_text_overall)
            )),
          generationConfig = list(
            temperature = temperature,
            maxOutputTokens = max_output_tokens,
            seed = 123456  # For reproducible runs
          )
        )
      )
      
      if(response_overall$status_code>200) {
        stop(paste("Error - ", content(response_overall)$error$message))
      }
      
      # Extract the Gemini model's response
      candidates_overall <- httr::content(response_overall)$candidates
      outputs_overall <- unlist(lapply(candidates_overall, function(candidate) candidate$content$parts))
      
      # Print the overall summary
      cat("\nOverall Summary:\n")
      cat(outputs_overall)
      
      return(vis_network_combined)
    } else if (do_kegg) {
      # Save the network data to a tab-delimited text file
      write.table(network_df_kegg, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create an interactive network visualization
      network_df_kegg <- data.frame(
        source = network_df_kegg$Gene,
        target = network_df_kegg$Concept,
        interaction = network_df_kegg$Interaction
      )
      
      # Create the graph object
      network_graph_kegg <- igraph::graph_from_data_frame(network_df_kegg, directed = TRUE)
      
      # Create the visNetwork object
      vis_network_kegg <- visNetwork(nodes = data.frame(id = igraph::V(network_graph_kegg)$name, 
                                                        label = igraph::V(network_graph_kegg)$name), 
                                     edges = data.frame(from = network_df_kegg$source, 
                                                        to = network_df_kegg$target, 
                                                        label = network_df_kegg$interaction))
      
      # Customize the visualization
      vis_network_kegg <- vis_network_kegg %>%
        visEdges(arrows = "to") %>% 
        visNodes(size = 20,  # Adjust node size as needed
                 font = list(size = 24)) %>%  # Increase node label size
        visOptions(highlightNearest = TRUE)  # Highlight nodes when hovering
      
      # Save the visualization as an HTML file
      visSave(vis_network_kegg, file = html_file)
      return(vis_network_kegg)
    } else if (do_go) {
      # Save the network data to a tab-delimited text file
      write.table(network_df_go, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create an interactive network visualization
      network_df_go <- data.frame(
        source = network_df_go$Gene,
        target = network_df_go$Concept,
        interaction = network_df_go$Interaction
      )
      
      # Create the graph object
      network_graph_go <- igraph::graph_from_data_frame(network_df_go, directed = TRUE)
      
      # Create the visNetwork object
      vis_network_go <- visNetwork(nodes = data.frame(id = igraph::V(network_graph_go)$name, 
                                                      label = igraph::V(network_graph_go)$name), 
                                   edges = data.frame(from = network_df_go$source, 
                                                      to = network_df_go$target, 
                                                      label = network_df_go$interaction))
      
      # Customize the visualization
      vis_network_go <- vis_network_go %>%
        visEdges(arrows = "to") %>% 
        visNodes(size = 20,  # Adjust node size as needed
                 font = list(size = 24)) %>%  # Increase node label size
        visOptions(highlightNearest = TRUE)  # Highlight nodes when hovering
      
      # Save the visualization as an HTML file
      visSave(vis_network_go, file = html_file)
      return(vis_network_go)
    }
  })
}

SCassist_analyze_enrichment_L <- function(organism = "human", markers = "", 
                                          pvalue = 0.05, log2FC = 1, 
                                          experimental_design = "Single-cell RNA sequencing",
                                          output_file = "L_SCassist_summary_network_data.txt", 
                                          html_file = "L_SCassist_summary_network_visualization.html",
                                          do_kegg = TRUE, do_go = TRUE, go_ont = "BP",
                                          model = "llama3",
                                          model_params = list(seed = 42, temperature = 0, num_gpu = 0)) {
  
  geneID <- NULL
  . <- NULL
  Description <- NULL
  Gene <- NULL
  SYMBOL <- NULL
  subcategory <- NULL
  
  # Suppress the warning message
  suppressWarnings({
    # Function to check and install the database if needed
    install_organism_db <- function(organism) {
      organism_db <- switch(organism,
                            "human" = "org.Hs.eg.db",
                            "mouse" = "org.Mm.eg.db",
                            "rat" = "org.Rn.eg.db",
                            "zebrafish" = "org.Dr.eg.db",
                            "fly" = "org.Dm.eg.db",
                            "yeast" = "org.Sc.sgd.db",
                            "worm" = "org.Ce.eg.db",
                            "arabidopsis" = "org.At.tair.db",
                            "chicken" = "org.Gg.eg.db",
                            stop("Unsupported organism. Please choose from the list.")
      )
      
      if (!requireNamespace(organism_db, quietly = TRUE)) {
        cat(paste("The", organism_db, "database is not installed. Installing...\n"))
        BiocManager::install(organism_db)
      }
      
      library(organism_db, character.only = TRUE)
      return(organism_db)
    }
    
    # Convert organism name to a standard code used by clusterProfiler
    organism_code <- switch(organism,
                            "human" = "hsa",
                            "mouse" = "mmu",
                            "rat" = "rno",
                            "zebrafish" = "dre",
                            "fly" = "dme",
                            "yeast" = "sce",
                            "worm" = "cel",
                            "arabidopsis" = "ath",
                            "chicken" = "gga",
                            stop("Unsupported organism. Please choose from the list.")
    )
    
    # Check and install the database if needed
    organism_db <- install_organism_db(organism)
    
    # 2. Use :: to access the functions from the specific packages
    markers <- markers %>%
      dplyr::filter(markers[["p_val"]] < pvalue & (abs(markers[["avg_log2FC"]]) > log2FC))
    
    # 3. Extract gene names as a list
    gene_list <- rownames(markers) 
    gene_list <- sub("\\..*", "", gene_list)
    
    # 4. Determine input type (SYMBOL or ENTREZID)
    input_type <- ifelse(is.character(gene_list[[1]]), "SYMBOL", "ENTREZID")
    
    # 5. Create a table with gene names and IDs based on input type
    if (input_type == "ENTREZID") {
      gene_mapping <- suppressMessages(clusterProfiler::bitr(gene_list, fromType = "ENTREZID", toType = c("SYMBOL", "ENTREZID"), OrgDb = organism_db))
    } else {
      gene_mapping <- suppressMessages(clusterProfiler::bitr(gene_list, fromType = "SYMBOL", toType = c("SYMBOL", "ENTREZID"), OrgDb = organism_db))
    }
    
    # KEGG Enrichment
    if (do_kegg) {
      
      writeLines("\nThe KEGG analysis is presented below.
-------------------------------------------\n")
      
      # 6. Perform pathway enrichment analysis using ENTREZID
      kegg_enrichment <- clusterProfiler::enrichKEGG(gene = gene_mapping$ENTREZID, organism = organism_code)
      
      # 7. Get the output of kegg enrichment analysis into a separate table
      kegg_enrichment_results <- kegg_enrichment@result
 
      # Save enrichment results to a text file
      write.table(kegg_enrichment_results, file = "L_kegg_enrichment_results_original.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The ClusterProfiler KEGG pathways result is saved as a text file in the current working directory. \n")
      
      # 8. Filter for pvalue < 0.05
      kegg_enrichment_results <- kegg_enrichment_results %>%
        dplyr::filter(pvalue < 0.05)
      
      # 9. Map gene names to gene IDs
      kegg_enrichment_results <- kegg_enrichment_results %>%
        dplyr::mutate(geneID = strsplit(geneID, "/")) %>% 
        tidyr::unnest(geneID) %>%
        # Convert geneID to character for joining
        dplyr::mutate(geneID = as.character(geneID)) %>%
        # Join with gene_mapping based on geneID and ENTREZID
        dplyr::left_join(gene_mapping, by = c("geneID" = "ENTREZID")) %>%
        # Replace geneID with SYMBOL if available
        dplyr::mutate(Gene = ifelse(!is.na(SYMBOL), SYMBOL, geneID)) %>% 
        # Group by subcategory, Description, and pvalue
        dplyr::group_by(subcategory, Description, pvalue) %>%
        # Summarize the results, keeping only unique combinations
        dplyr::summarize(., geneID = paste(unique(geneID), collapse = ", "),
                         Gene = paste(unique(Gene), collapse = ", "), .groups = "drop") 
      
      # Alternatively, use the square bracket notation
      input_kegg_enrichment <- kegg_enrichment_results[, c("Description", "pvalue", "Gene")]
      
      # 13. Save enrichment results to a text file
      write.table(kegg_enrichment_results, file = "L_kegg_enrichment_results_genename.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The list of significantly enriched KEGG pathways is also saved as a text file in the current working directory. \n")
      
      #json_kegg_enrichment <- jsonlite::toJSON(input_kegg_enrichment)
      
      input_text_kegg <- paste0("
This below Data is a list of KEGG pathway enrichment results for a set of differentially expressed genes obtained from a single cell experiment involving ",experimental_design,".\n 
Data:
```json \n",
                                input_kegg_enrichment,
                                "```\n
Analyze all of the pathways from my above Data and provide insights in a structured format. Include:
                           
        1. **Significant Pathways:** Analyze all the pathways in my list, in the context of the system involving ",experimental_design," and summarize common themes, in a small paragraph.
        2. **Regulators:** Include potential involvment of any transcription factors from the given list, in the context of the system involving ",experimental_design," in a 5 line paragraph.
        3. **Key Genes or Potential Targets:** Suggest which genes from the enriched pathways, in my data input, might be important to the system or a target based on their potential impact on the system involving", experimental_design," in a 10 line paragraph.
        
        Do not provide any Further Investigation or comments, not asked for.")
      
      response_kegg <- suppressMessages(rollama::query(
        input_text_kegg,
        screen = FALSE,
        model = model,
        model_params = model_params))
      
      cat(response_kegg[[1]]$message$content)
      
      # 14. Return the raw output from the Gemini Pro language model
      input_for_network_kegg <- paste0("Extract important named entities and relationships between key genes and 
                                potential concepts, not cell types, and format it in three columns as Gene, 
                                Interaction, Concept, one gene - one interaction per row. 
                                the gene column should only contain one gene. Include only concepts related to biological pathways.
                                the concept column should only contain ONE gene or ONE concept. 
                                Add potential gene-gene interactions as well, for the core system. 
                                Do not provide any concept - concept associations.
                                EXAMPLE OUTPUT:
                                | Gene1 | Involved in | Metabolism |
                                | Gene3 | Interacts with | Gene5 |\n",
                                       response_kegg[[1]]$message$content)
      
      response2_kegg <- suppressMessages(rollama::query(
        input_for_network_kegg,
        screen = FALSE,
        model = model,
        model_params = model_params))
      
      # Extract and format the network data from the LLM response
      network_section_kegg <- strsplit(response2_kegg[[1]]$message$content, "\n")[[1]]
      network_data_kegg <- network_section_kegg[5:(length(network_section_kegg)-2)]
      network_data_kegg <- strsplit(network_data_kegg, "\\|\\s*")
      
      # Create a data frame from the network data
      network_df_kegg <- data.frame(
        Gene = sapply(network_data_kegg, "[", 2),
        Interaction = sapply(network_data_kegg, "[", 3),
        Concept = sapply(network_data_kegg, "[", 4)
      )
      
      # Clean up the data frame
      network_df_kegg <- data.frame(
        lapply(network_df_kegg, function(x) gsub("\\s+", " ", trimws(x)))
      )
    }
    
    # GO Enrichment
    if (do_go) {
      
      writeLines("\n\nThe GO analysis is provided below.
----------------------------------------\n")
      
      # 6. Perform pathway enrichment analysis using ENTREZID
      go_enrichment <- clusterProfiler::enrichGO(gene = gene_mapping$ENTREZID, OrgDb = organism_db, ont=go_ont)
      
      # 7. Get the output of kegg enrichment analysis into a separate table
      go_enrichment_results <- go_enrichment@result

      # Save enrichment results to a text file
      write.table(go_enrichment_results, file = "L_go_enrichment_results_original.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The ClusterProfiler GO enrichment result is saved as a text file in the current working directory. \n")
      
      # 8. Filter for pvalue < 0.05
      go_enrichment_results <- go_enrichment_results %>%
        dplyr::filter(pvalue < 0.05)
      
      # 9. Map gene names to gene IDs
      go_enrichment_results <- go_enrichment_results %>%
        dplyr::mutate(geneID = strsplit(geneID, "/")) %>% 
        tidyr::unnest(geneID) %>%
        # Convert geneID to character for joining
        dplyr::mutate(geneID = as.character(geneID)) %>%
        # Join with gene_mapping based on geneID and ENTREZID
        dplyr::left_join(gene_mapping, by = c("geneID" = "ENTREZID")) %>%
        # Replace geneID with SYMBOL if available
        dplyr::mutate(Gene = ifelse(!is.na(SYMBOL), SYMBOL, geneID)) %>% 
        # Group by subcategory, Description, and pvalue
        dplyr::group_by(Description, pvalue) %>%
        # Summarize the results, keeping only unique combinations
        dplyr::summarize(., geneID = paste(unique(geneID), collapse = ", "),
                         Gene = paste(unique(Gene), collapse = ", "), .groups = "drop") 
      
      # Alternatively, use the square bracket notation
      input_go_enrichment <- go_enrichment_results[, c("Description", "pvalue", "Gene")]
 
      # Save enrichment results to a text file
      write.table(go_enrichment_results, file = "L_go_enrichment_results_genename.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The list of significantly enriched GO terms is also saved as a text file in the current working directory. \n")
      
      #json_go_enrichment <- jsonlite::toJSON(input_go_enrichment)
      
      input_text_go <- paste0("
This below Data is a list of GO enrichment results for a set of differentially expressed genes obtained from a single cell experiment involving ",experimental_design,".\n 
Data:
```json \n",
                              input_go_enrichment,
                              "```\n
Analyze all of the enriched concepts from my above Data and provide insights in a structured format. Include:
                           
        1. **Significant Concepts:** Analyze all the concepts in my list, in the context of the system involving ",experimental_design," and summarize common themes, in 2 paragraph, with a total of no more than 10 lines.
        
        2. **Regulators:** Include potential involvment of any transcription factors from the given list, in the context of the system involving ",experimental_design," in a 5 line paragraph.
        
        3. **Key Genes or Potential Targets:** Suggest which genes from the enriched ontologies, in my data input, might be important to the system or a target based on their potential impact on the system involving", experimental_design," in a 10 line paragraph.
        
        Do not provide any Further Investigation or comments, not asked for.")
      
      response_go <- suppressMessages(rollama::query(
        input_text_go,
        screen = FALSE,
        model = model,
        model_params = model_params))
      
      cat(response_go[[1]]$message$content)
      
      # 14. Return the raw output from the Gemini Pro language model
      input_for_network_go <- paste0("Extract important named entities and relationships between key genes and 
                                potential concepts, not cell types, and format it in three columns as Gene, 
                                Interaction, Concept, one gene - one interaction per row. 
                                the gene column should only contain one gene. Include only concepts related to GO terms.
                                the concept column should only contain ONE gene or ONE concept. 
                                Add potential gene-gene interactions as well, for the core system. 
                                Do not provide any concept - concept associations.
                                EXAMPLE OUTPUT:
                                | Gene1 | Involved in | Metabolism |
                                | Gene3 | Interacts with | Gene5 |\n",
                                     response_go[[1]]$message$content)
      
      response2_go <- suppressMessages(rollama::query(
        input_for_network_go,
        screen = FALSE,
        model = model,
        model_params = model_params))
      
      # Extract and format the network data from the LLM response
      network_section_go <- strsplit(response2_go[[1]]$message$content, "\n")[[1]]
      network_data_go <- network_section_go[5:(length(network_section_go)-2)]
      network_data_go <- strsplit(network_data_go, "\\|\\s*")
      
      # Create a data frame from the network data
      network_df_go <- data.frame(
        Gene = sapply(network_data_go, "[", 2),
        Interaction = sapply(network_data_go, "[", 3),
        Concept = sapply(network_data_go, "[", 4)
      )
      
      # Clean up the data frame
      network_df_go <- data.frame(
        lapply(network_df_go, function(x) gsub("\\s+", " ", trimws(x)))
      )
    }
    
    # Combine Networks
    if (do_kegg && do_go) {
      network_df_combined <- rbind(network_df_kegg, network_df_go)
      
      # Save the combined network data to a tab-delimited text file
      write.table(network_df_combined, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create an interactive network visualization
      network_df_combined <- data.frame(
        source = network_df_combined$Gene,
        target = network_df_combined$Concept,
        interaction = network_df_combined$Interaction
      )
      
      # Create the graph object
      network_graph_combined <- igraph::graph_from_data_frame(network_df_combined, directed = TRUE)
      
      # Create the visNetwork object
      vis_network_combined <- visNetwork(nodes = data.frame(id = igraph::V(network_graph_combined)$name, 
                                                            label = igraph::V(network_graph_combined)$name), 
                                         edges = data.frame(from = network_df_combined$source, 
                                                            to = network_df_combined$target, 
                                                            label = network_df_combined$interaction))
      
      # Customize the visualization
      vis_network_combined <- vis_network_combined %>%
        visEdges(arrows = "to") %>% 
        visNodes(size = 20,  # Adjust node size as needed
                 font = list(size = 24)) %>%  # Increase node label size
        visOptions(highlightNearest = TRUE)  # Highlight nodes when hovering
      
      # Save the visualization as an HTML file
      visSave(vis_network_combined, file = html_file)
      
      # Combine Gemini Pro responses for overall summary
      combined_response <- paste0(response_kegg$message$content, "\n\n", response_go$message$content)
      
      # Generate overall summary prompt
      input_text_overall <- paste0("
This is a combined summary of KEGG and GO enrichment results for a set of differentially expressed genes obtained from a single cell experiment involving ",experimental_design,".\n 
Please provide a comprehensive summary based on the following information:
", combined_response, "\n
Do not provide any Further Investigation or comments, not asked for.")
      
      response_overall <- suppressMessages(rollama::query(
        input_text_overall,
        screen = FALSE,
        model = model,
        model_params = model_params))
      
      # Print the overall summary
      cat("\nOverall Summary:\n")
      cat(response_overall[[1]]$message$content)
      
      return(vis_network_combined)
    } else if (do_kegg) {
      # Save the network data to a tab-delimited text file
      write.table(network_df_kegg, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create an interactive network visualization
      network_df_kegg <- data.frame(
        source = network_df_kegg$Gene,
        target = network_df_kegg$Concept,
        interaction = network_df_kegg$Interaction
      )
      
      # Create the graph object
      network_graph_kegg <- igraph::graph_from_data_frame(network_df_kegg, directed = TRUE)
      
      # Create the visNetwork object
      vis_network_kegg <- visNetwork(nodes = data.frame(id = igraph::V(network_graph_kegg)$name, 
                                                        label = igraph::V(network_graph_kegg)$name), 
                                     edges = data.frame(from = network_df_kegg$source, 
                                                        to = network_df_kegg$target, 
                                                        label = network_df_kegg$interaction))
      
      # Customize the visualization
      vis_network_kegg <- vis_network_kegg %>%
        visEdges(arrows = "to") %>% 
        visNodes(size = 20,  # Adjust node size as needed
                 font = list(size = 24)) %>%  # Increase node label size
        visOptions(highlightNearest = TRUE)  # Highlight nodes when hovering
      
      # Save the visualization as an HTML file
      visSave(vis_network_kegg, file = html_file)
      return(vis_network_kegg)
    } else if (do_go) {
      # Save the network data to a tab-delimited text file
      write.table(network_df_go, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create an interactive network visualization
      network_df_go <- data.frame(
        source = network_df_go$Gene,
        target = network_df_go$Concept,
        interaction = network_df_go$Interaction
      )
      
      # Create the graph object
      network_graph_go <- igraph::graph_from_data_frame(network_df_go, directed = TRUE)
      
      # Create the visNetwork object
      vis_network_go <- visNetwork(nodes = data.frame(id = igraph::V(network_graph_go)$name, 
                                                      label = igraph::V(network_graph_go)$name), 
                                   edges = data.frame(from = network_df_go$source, 
                                                      to = network_df_go$target, 
                                                      label = network_df_go$interaction))
      
      # Customize the visualization
      vis_network_go <- vis_network_go %>%
        visEdges(arrows = "to") %>% 
        visNodes(size = 20,  # Adjust node size as needed
                 font = list(size = 24)) %>%  # Increase node label size
        visOptions(highlightNearest = TRUE)  # Highlight nodes when hovering
      
      # Save the visualization as an HTML file
      visSave(vis_network_go, file = html_file)
      return(vis_network_go)
    }
  })
}
















SCassist_analyze_enrichment_C <- function(organism = "human", markers="", 
                                          pvalue = 0.05, log2FC = 1, 
                                          experimental_design = "Single-cell RNA sequencing",
                                          temperature = 0,
                                          max_output_tokens = 10048,
                                          model = "gpt-4o-mini",
                                          model_params = list(seed = 42, temperature = 0),
                                          api_key_file = api_key_file,
                                          output_file = "C_SCassist_summary_network_data.txt", 
                                          html_file = "C_SCassist_summary_network_visualization.html",
                                          do_kegg = TRUE, do_go = TRUE, go_ont = "BP") {
  
  Description <- Gene <- SYMBOL <- geneID <- V <- subcategory <- . <- NULL
  
  # 1. Read the API key from the specified file
  api_key <- readLines(api_key_file, encoding = "UTF-8")
  
  # Suppress the warning message
  suppressWarnings({
    # Function to check and install the database if needed
    install_organism_db <- function(organism) {
      organism_db <- switch(organism,
                            "human" = "org.Hs.eg.db",
                            "mouse" = "org.Mm.eg.db",
                            "rat" = "org.Rn.eg.db",
                            "zebrafish" = "org.Dr.eg.db",
                            "fly" = "org.Dm.eg.db",
                            "yeast" = "org.Sc.sgd.db",
                            "worm" = "org.Ce.eg.db",
                            "arabidopsis" = "org.At.tair.db",
                            "chicken" = "org.Gg.eg.db",
                            stop("Unsupported organism. Please choose from the list.")
      )
      
      if (!requireNamespace(organism_db, quietly = TRUE)) {
        cat(paste("The", organism_db, "database is not installed. Installing...\n"))
        BiocManager::install(organism_db)
      }
      
      # Use :: to access the library function from the specific package
      library(organism_db, character.only = TRUE)
      return(organism_db)
    }
    
    # Convert organism name to a standard code used by clusterProfiler
    organism_code <- switch(organism,
                            "human" = "hsa",
                            "mouse" = "mmu",
                            "rat" = "rno",
                            "zebrafish" = "dre",
                            "fly" = "dme",
                            "yeast" = "sce",
                            "worm" = "cel",
                            "arabidopsis" = "ath",
                            "chicken" = "gga",
                            stop("Unsupported organism. Please choose from the list.")
    )
    
    # 1. Check and install the database if needed
    organism_db <- install_organism_db(organism)
    
    # 2. Use :: to access the functions from the specific packages
    markers <- markers %>%
      dplyr::filter(markers[["p_val"]] < pvalue & (abs(markers[["avg_log2FC"]]) > log2FC))
    
    # 3. Extract gene names as a list
    gene_list <- rownames(markers) 
    gene_list <- sub("\\..*", "", gene_list)
    
    # 4. Determine input type (SYMBOL or ENTREZID)
    input_type <- ifelse(is.character(gene_list[[1]]), "SYMBOL", "ENTREZID")
    
    # 5. Create a table with gene names and IDs based on input type
    if (input_type == "ENTREZID") {
      gene_mapping <- suppressMessages(clusterProfiler::bitr(gene_list, fromType = "ENTREZID", toType = c("SYMBOL", "ENTREZID"), OrgDb = organism_db))
    } else {
      gene_mapping <- suppressMessages(clusterProfiler::bitr(gene_list, fromType = "SYMBOL", toType = c("SYMBOL", "ENTREZID"), OrgDb = organism_db))
    }
    
    # KEGG Enrichment
    if (do_kegg) {
      
      writeLines("\nThe KEGG analysis is presented below.
-------------------------------------------\n")
      
      # 6. Perform pathway enrichment analysis using ENTREZID
      kegg_enrichment <- clusterProfiler::enrichKEGG(gene = gene_mapping$ENTREZID, organism = organism_code)
      
      # 7. Get the output of kegg enrichment analysis into a separate table
      kegg_enrichment_results <- kegg_enrichment@result
      
      # 13. Save enrichment results to a text file
      write.table(kegg_enrichment_results, file = "C_kegg_enrichment_results_original.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The ClusterProfiler KEGG pathways result is saved as a text file in the current working directory. \n")
      
      # 8. Filter for pvalue < 0.05
      kegg_enrichment_results <- kegg_enrichment_results %>%
        dplyr::filter(pvalue < 0.05)
      
      # 9. Map gene names to gene IDs
      kegg_enrichment_results <- kegg_enrichment_results %>%
        dplyr::mutate(geneID = strsplit(geneID, "/")) %>% 
        tidyr::unnest(geneID) %>%
        # Convert geneID to character for joining
        dplyr::mutate(geneID = as.character(geneID)) %>%
        # Join with gene_mapping based on geneID and ENTREZID
        dplyr::left_join(gene_mapping, by = c("geneID" = "ENTREZID")) %>%
        # Replace geneID with SYMBOL if available
        dplyr::mutate(Gene = ifelse(!is.na(SYMBOL), SYMBOL, geneID)) %>% 
        # Group by subcategory, Description, and pvalue
        dplyr::group_by(subcategory, Description, pvalue) %>%
        # Summarize the results, keeping only unique combinations
        dplyr::summarize(., geneID = paste(unique(geneID), collapse = ", "),
                         Gene = paste(unique(Gene), collapse = ", "), .groups = "drop") 
      
      # Use the square bracket notation
      input_kegg_enrichment <- kegg_enrichment_results[, c("Description", "pvalue", "Gene")]
      
      # 13. Save enrichment results to a text file
      write.table(kegg_enrichment_results, file = "C_kegg_enrichment_results_genename.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The list of significantly enriched KEGG pathways is also saved as a text file in the current working directory. \n")
      
      json_kegg_enrichment <- jsonlite::toJSON(input_kegg_enrichment)
      
      input_text_kegg <- paste0("
This below Data is a list of KEGG pathway enrichment results for a set of differentially expressed genes obtained from a single cell experiment involving ",experimental_design,".\n 
Data:
```json \n",
                                json_kegg_enrichment,
                                "```\n
Analyze all of the pathways from my above Data and provide insights in a structured format. Include:
                           
        1. **Significant Pathways:** Analyze all the pathways in my list, in the context of the system involving ",experimental_design," and summarize common themes, in 2 paragraph, with a total of 10 lines.
        
        2. **Regulators:** Include potential involvment of any transcription factors from the given list, in the context of the system involving ",experimental_design," in a 5 line paragraph.
        
        3. **Key Genes or Potential Targets:** Suggest which genes from the enriched pathways, in my data input, might be important to the system or a target based on their potential impact on the system involving", experimental_design," in a 10 line paragraph.
        
        Do not provide any Further Investigation or comments, not asked for.")
      
      # Make the POST request to OpenAI
      response_kegg <- httr::POST(
        url = "https://api.openai.com/v1/chat/completions",
        add_headers("Authorization" = paste("Bearer", api_key)),
        content_type_json(), 
        encode = "json",
        body = list(
          model = model, # Specify the model
          messages = list(      # Messages array as required by Chat Completions API
            list(
              role = "user",
              content = input_text_kegg
            )
          )
        )
      )
      
      # 12. Extract the OpenAI model's response
      # Check for HTTP errors
      if (http_error(response_kegg)) {
        stop(paste("OpenAI API request failed with status", http_status(response_kegg)$message, "\nContent:", content(response_kegg, "text")))
      }
      
      # Parse the JSON response
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required but not installed.")
      }
      
      response_content <- httr::content(response_kegg, "text")
      response_json <- jsonlite::fromJSON(response_content, flatten = TRUE)
      
      
      # Extract and print the generated text from the response
      # For chat completions, the text is in choices[[1]]$message$content
      if (!is.null(response_json$choices) && length(response_json$choices) > 0 && !is.null(response_json$choices$message.content)) {
        generated_text <- response_json$choices$message.content
        #cat(generated_text)
      } else {
        cat("Error: Could not extract generated text from OpenAI response.\n")
        print(response_json) # Print the full response for debugging
      }
      
      # 13. Print and return the LLM's response
      #cat(generated_text)
      outputs_kegg <- generated_text
      
      # 12. Print the LLM's response
      cat("KEGG Enrichment Insights:\n")
      cat(outputs_kegg)
      
      # 14. Return the raw output from the Gemini Pro language model
      input_for_network_kegg <- paste0("Extract important named entities and relationships between key genes and 
                                potential concepts, not cell types, and format it in three columns as Gene, 
                                Interaction, Concept, one gene - one interaction per row. 
                                the gene column should only contain one gene. Include only concepts related to biological pathways.
                                the concept column should only contain ONE gene or ONE concept. 
                                Add potential gene-gene interactions as well, for the core system. 
                                Do not provide any concept - concept associations.
                                EXAMPLE OUTPUT:
                                | Gene1 | Involved in | Metabolism |
                                | Gene3 | Interacts with | Gene5 |\n",
                                       outputs_kegg)
      
      # Make the POST request to OpenAI
      response_kegg_network <- httr::POST(
        url = "https://api.openai.com/v1/chat/completions",
        add_headers("Authorization" = paste("Bearer", api_key)),
        content_type_json(), 
        encode = "json",
        body = list(
          model = model, # Specify the model
          messages = list(      # Messages array as required by Chat Completions API
            list(
              role = "user",
              content = input_for_network_kegg
            )
          )
        )
      )
      # 12. Extract the OpenAI model's response
      # Check for HTTP errors
      if (http_error(response_kegg_network)) {
        stop(paste("OpenAI API request failed with status", http_status(response_kegg_network)$message, "\nContent:", content(response_kegg_network, "text")))
      }
      
      # Parse the JSON response
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required but not installed.")
      }
      
      response_content <- httr::content(response_kegg_network, "text")
      response_json <- jsonlite::fromJSON(response_content, flatten = TRUE)
      
      
      # Extract and print the generated text from the response
      # For chat completions, the text is in choices[[1]]$message$content
      if (!is.null(response_json$choices) && length(response_json$choices) > 0 && !is.null(response_json$choices$message.content)) {
        generated_text <- response_json$choices$message.content
        #cat(generated_text)
      } else {
        cat("Error: Could not extract generated text from OpenAI response.\n")
        print(response_json) # Print the full response for debugging
      }
      
      # 13. Print and return the LLM's response
      #cat(generated_text)
      outputs2_kegg <- generated_text
      
      # 12. Print the LLM's response
      
      # Extract and format the network data from the LLM response
      network_section_kegg <- strsplit(outputs2_kegg, "\n")[[1]]
      network_data_kegg <- network_section_kegg[5:length(network_section_kegg)]
      network_data_kegg <- strsplit(network_data_kegg, "\\|\\s*")
      
      #cat(paste(network_data_kegg, collapse = "\n"))
      #str(network_data_kegg) 
      #class(network_data_kegg)
      
      # Create a data frame from the network data
      network_df_kegg <- data.frame(
        Gene = sapply(network_data_kegg, "[", 2),
        Interaction = sapply(network_data_kegg, "[", 3),
        Concept = sapply(network_data_kegg, "[", 4)
      )
      
      # Clean up the data frame
      network_df_kegg <- data.frame(
        lapply(network_df_kegg, function(x) gsub("\\s+", " ", trimws(x)))
      )
    }
    
    # GO Enrichment
    if (do_go) {
      
      writeLines("\nThe GO analysis is provided below.
----------------------------------------\n")
      
      # 6. Perform pathway enrichment analysis using ENTREZID
      go_enrichment <- clusterProfiler::enrichGO(gene = gene_mapping$ENTREZID, OrgDb = organism_db, ont=go_ont)
      
      # 7. Get the output of kegg enrichment analysis into a separate table
      go_enrichment_results <- go_enrichment@result
      
      # 13. Save enrichment results to a text file
      write.table(go_enrichment_results, file = "C_go_enrichment_results_original.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The ClusterProfiler GO enrichment result is saved as a text file in the current working directory. \n")
      
      # 8. Filter for pvalue < 0.05, keep only top 100 for openai context limit
      go_enrichment_results <- go_enrichment_results %>%
        dplyr::filter(pvalue < 0.005) %>%
        dplyr::arrange(pvalue) %>%
        head(200)
      
      # 9. Map gene names to gene IDs
      go_enrichment_results <- go_enrichment_results %>%
        dplyr::mutate(geneID = strsplit(geneID, "/")) %>% 
        tidyr::unnest(geneID) %>%
        # Convert geneID to character for joining
        dplyr::mutate(geneID = as.character(geneID)) %>%
        # Join with gene_mapping based on geneID and ENTREZID
        dplyr::left_join(gene_mapping, by = c("geneID" = "ENTREZID")) %>%
        # Replace geneID with SYMBOL if available
        dplyr::mutate(Gene = ifelse(!is.na(SYMBOL), SYMBOL, geneID)) %>% 
        # Group by subcategory, Description, and pvalue
        dplyr::group_by(Description, pvalue) %>%
        # Summarize the results, keeping only unique combinations
        dplyr::summarize(., geneID = paste(unique(geneID), collapse = ", "),
                         Gene = paste(unique(Gene), collapse = ", "), .groups = "drop") 
      
      # Alternatively, use the square bracket notation
      input_go_enrichment <- go_enrichment_results[, c("Description", "pvalue", "Gene")]
      
      # 13. Save enrichment results to a text file
      write.table(go_enrichment_results, file = "C_go_enrichment_results_genename.txt", sep = "\t",
                  row.names = FALSE, quote = FALSE)
      
      writeLines("The list of significantly enriched GO terms is also saved as a text file in the current working directory. \n")
      
      json_go_enrichment <- jsonlite::toJSON(input_go_enrichment)
      
      input_text_go <- paste0("
This below Data is a list of GO enrichment results for a set of differentially expressed genes obtained from a single cell experiment involving ",experimental_design,".\n 
Data:
```json \n",
                              json_go_enrichment,
                              "```\n
Analyze all of the enriched concepts from my above Data and provide insights in a structured format. Include:
                           
        1. **Significant Concepts:** Analyze all the concepts in my list, in the context of the system involving ",experimental_design," and summarize common themes, in 2 paragraph, with a total of no more than 10 lines.
        
        2. **Regulators:** Include potential involvment of any transcription factors from the given list, in the context of the system involving ",experimental_design," in a 5 line paragraph. Restrict the number of genes to less than 10.
        
        3. **Key Genes or Potential Targets:** Suggest which genes from the enriched ontologies, in my data input, might be important to the system or a target based on their potential impact on the system involving", experimental_design," in a 10 line paragraph.
        
        Do not provide any Further Investigation or comments, not asked for.")
      
      # Make the POST request to OpenAI
      response <- httr::POST(
        url = "https://api.openai.com/v1/chat/completions",
        add_headers("Authorization" = paste("Bearer", api_key)),
        content_type_json(), 
        encode = "json",
        body = list(
          model = model, # Specify the model
          messages = list(      # Messages array as required by Chat Completions API
            list(
              role = "user",
              content = input_text_go
            )
          )
        )
      )
      
      # 12. Extract the OpenAI model's response
      # Check for HTTP errors
      if (http_error(response)) {
        stop(paste("OpenAI API request failed with status", http_status(response)$message, "\nContent:", httr::content(response, "text")))
      }
      
      # Parse the JSON response
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required but not installed.")
      }
      
      response_content <- httr::content(response, "text")
      response_json <- jsonlite::fromJSON(response_content, flatten = TRUE)
      
      
      # Extract and print the generated text from the response
      # For chat completions, the text is in choices[[1]]$message$content
      if (!is.null(response_json$choices) && length(response_json$choices) > 0 && !is.null(response_json$choices$message.content)) {
        generated_text <- response_json$choices$message.content
        #cat(generated_text)
      } else {
        cat("Error: Could not extract generated text from OpenAI response.\n")
        print(response_json) # Print the full response for debugging
      }
      
      # 13. Print and return the LLM's response
      #cat(generated_text)
      outputs_go <- generated_text
      
      # 12. Print the LLM's response
      cat("GO Enrichment Insights:\n")
      cat(outputs_go)
      
      # 14. Return the raw output from the Gemini Pro language model
      input_for_network_go <- paste0("Extract important named entities and relationships between key genes and 
                                potential concepts, not cell types, and format it in three columns as Gene, 
                                Interaction, Concept, one gene - one interaction per row. 
                                the gene column should only contain one gene. Include only concepts related to GO terms.
                                the concept column should only contain ONE gene or ONE concept. 
                                Add potential gene-gene interactions as well, for the core system. 
                                Do not provide any concept - concept associations.
                                EXAMPLE OUTPUT:
                                | Gene1 | Involved in | Metabolism |
                                | Gene3 | Interacts with | Gene5 |\n",
                                     outputs_go)
      
      # Make the POST request to OpenAI
      response <- httr::POST(
        url = "https://api.openai.com/v1/chat/completions",
        add_headers("Authorization" = paste("Bearer", api_key)),
        content_type_json(), 
        encode = "json",
        body = list(
          model = model, # Specify the model
          messages = list(      # Messages array as required by Chat Completions API
            list(
              role = "user",
              content = input_for_network_go
            )
          )
        )
      )
      
      # 12. Extract the OpenAI model's response
      # Check for HTTP errors
      if (http_error(response)) {
        stop(paste("OpenAI API request failed with status", http_status(response)$message, "\nContent:", httr::content(response, "text")))
      }
      
      # Parse the JSON response
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required but not installed.")
      }
      
      response_content <- httr::content(response, "text")
      response_json <- jsonlite::fromJSON(response_content, flatten = TRUE)
      
      
      # Extract and print the generated text from the response
      # For chat completions, the text is in choices[[1]]$message$content
      if (!is.null(response_json$choices) && length(response_json$choices) > 0 && !is.null(response_json$choices$message.content)) {
        generated_text <- response_json$choices$message.content
        #cat(generated_text)
      } else {
        cat("Error: Could not extract generated text from OpenAI response.\n")
        print(response_json) # Print the full response for debugging
      }
      
      # 13. Print and return the LLM's response
      #cat(generated_text)
      
      outputs2_go <- generated_text
      
      # 12. Print the LLM's response
      
      # Extract and format the network data from the LLM response
      network_section_go <- strsplit(outputs2_go, "\n")[[1]]
      network_data_go <- network_section_go[5:length(network_section_go)]
      network_data_go <- strsplit(network_data_go, "\\|\\s*")
      
      # Create a data frame from the network data
      network_df_go <- data.frame(
        Gene = sapply(network_data_go, "[", 2),
        Interaction = sapply(network_data_go, "[", 3),
        Concept = sapply(network_data_go, "[", 4)
      )
      
      # Clean up the data frame
      network_df_go <- data.frame(
        lapply(network_df_go, function(x) gsub("\\s+", " ", trimws(x)))
      )
    }
    
    # Combine Networks
    if (do_kegg && do_go) {
      network_df_combined <- rbind(network_df_kegg, network_df_go)
      
      # Save the combined network data to a tab-delimited text file
      write.table(network_df_combined, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create an interactive network visualization
      network_df_combined <- data.frame(
        source = network_df_combined$Gene,
        target = network_df_combined$Concept,
        interaction = network_df_combined$Interaction
      )
      
      # Create the graph object
      network_graph_combined <- igraph::graph_from_data_frame(network_df_combined, directed = TRUE)
      
      # Create the visNetwork object
      vis_network_combined <- visNetwork(nodes = data.frame(id = igraph::V(network_graph_combined)$name, 
                                                            label = igraph::V(network_graph_combined)$name), 
                                         edges = data.frame(from = network_df_combined$source, 
                                                            to = network_df_combined$target, 
                                                            label = network_df_combined$interaction))
      
      # Customize the visualization
      vis_network_combined <- vis_network_combined %>%
        visEdges(arrows = "to") %>% 
        visNodes(size = 20,  # Adjust node size as needed
                 font = list(size = 24)) %>%  # Increase node label size
        visOptions(highlightNearest = TRUE)  # Highlight nodes when hovering
      
      # Save the visualization as an HTML file
      visSave(vis_network_combined, file = html_file)
      
      # Combine Gemini Pro responses for overall summary
      combined_response <- paste0(outputs_kegg, "\n\n", outputs_go)
      
      # Generate overall summary prompt
      input_text_overall <- paste0("
This is a combined summary of KEGG and GO enrichment results for a set of differentially expressed genes obtained from a single cell experiment involving ",experimental_design,".\n 
Please provide a comprehensive summary based on the following information:
", combined_response, "\n
Do not provide any Further Investigation or comments, not asked for.")
      
      # Make the POST request to OpenAI
      response <- httr::POST(
        url = "https://api.openai.com/v1/chat/completions",
        add_headers("Authorization" = paste("Bearer", api_key)),
        content_type_json(), 
        encode = "json",
        body = list(
          model = model, # Specify the model
          messages = list(      # Messages array as required by Chat Completions API
            list(
              role = "user",
              content = input_text_overall
            )
          )
        )
      )
      
      # 12. Extract the OpenAI model's response
      # Check for HTTP errors
      if (http_error(response)) {
        stop(paste("OpenAI API request failed with status", http_status(response)$message, "\nContent:", httr::content(response, "text")))
      }
      
      # Parse the JSON response
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required but not installed.")
      }
      
      response_content <- httr::content(response, "text")
      response_json <- jsonlite::fromJSON(response_content, flatten = TRUE)
      
      
      # Extract and print the generated text from the response
      # For chat completions, the text is in choices[[1]]$message$content
      if (!is.null(response_json$choices) && length(response_json$choices) > 0 && !is.null(response_json$choices$message.content)) {
        generated_text <- response_json$choices$message.content
        #cat(generated_text)
      } else {
        cat("Error: Could not extract generated text from OpenAI response.\n")
        print(response_json) # Print the full response for debugging
      }
      
      # 13. Print and return the LLM's response
      #cat(generated_text)

      outputs_overall <- generated_text
      
      # Print the overall summary
      cat("\nOverall Summary:\n")
      cat(outputs_overall)
      
      return(vis_network_combined)
    } else if (do_kegg) {
      # Save the network data to a tab-delimited text file
      write.table(network_df_kegg, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create an interactive network visualization
      network_df_kegg <- data.frame(
        source = network_df_kegg$Gene,
        target = network_df_kegg$Concept,
        interaction = network_df_kegg$Interaction
      )
      
      # Create the graph object
      network_graph_kegg <- igraph::graph_from_data_frame(network_df_kegg, directed = TRUE)
      
      # Create the visNetwork object
      vis_network_kegg <- visNetwork(nodes = data.frame(id = igraph::V(network_graph_kegg)$name, 
                                                        label = igraph::V(network_graph_kegg)$name), 
                                     edges = data.frame(from = network_df_kegg$source, 
                                                        to = network_df_kegg$target, 
                                                        label = network_df_kegg$interaction))
      
      # Customize the visualization
      vis_network_kegg <- vis_network_kegg %>%
        visEdges(arrows = "to") %>% 
        visNodes(size = 20,  # Adjust node size as needed
                 font = list(size = 24)) %>%  # Increase node label size
        visOptions(highlightNearest = TRUE)  # Highlight nodes when hovering
      
      # Save the visualization as an HTML file
      visSave(vis_network_kegg, file = html_file)
      return(vis_network_kegg)
    } else if (do_go) {
      # Save the network data to a tab-delimited text file
      write.table(network_df_go, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create an interactive network visualization
      network_df_go <- data.frame(
        source = network_df_go$Gene,
        target = network_df_go$Concept,
        interaction = network_df_go$Interaction
      )
      
      # Create the graph object
      network_graph_go <- igraph::graph_from_data_frame(network_df_go, directed = TRUE)
      
      # Create the visNetwork object
      vis_network_go <- visNetwork(nodes = data.frame(id = igraph::V(network_graph_go)$name, 
                                                      label = igraph::V(network_graph_go)$name), 
                                   edges = data.frame(from = network_df_go$source, 
                                                      to = network_df_go$target, 
                                                      label = network_df_go$interaction))
      
      # Customize the visualization
      vis_network_go <- vis_network_go %>%
        visEdges(arrows = "to") %>% 
        visNodes(size = 20,  # Adjust node size as needed
                 font = list(size = 24)) %>%  # Increase node label size
        visOptions(highlightNearest = TRUE)  # Highlight nodes when hovering
      
      # Save the visualization as an HTML file
      visSave(vis_network_go, file = html_file)
      return(vis_network_go)
    }
  })
}
