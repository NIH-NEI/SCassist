#' @title Analyze Top Principal Components (PCs) with a Large Language Model (LLM)
#'
#' @description This function analyzes the top principal components (PCs) from a Seurat object
#' using a large language model (LLM) to identify potential biological processes
#' driving the variations captured by each PC. The LLM considers the top genes
#' contributing to each PC and their associated pathways to provide insights
#' into the underlying biological system. This version analyzes a specified number of PCs
#' and then combines the summaries from the LLM for each PC into a single overall summary.
#'
#' @author Vijay Nagarajan, PhD, NEI/NIH
#'
#' @details This function was written with assistance from Google's Gemini.
#'
#' @param llm_server The LLM server to use. Options are "google" or "ollama" or "openai". Default is "google".
#' @param model_params A list of parameters to be passed to the `ollama::query` function.
#'   This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.
#' @param seurat_object_name The name of the Seurat object containing the
#'  single-cell RNA-seq data. The object should be accessible in the current
#'  environment and should have PCA already run (e.g., using `RunPCA`).
#' @param num_pcs The number of top PCs to analyze. Defaults to 6.
#' @param top_n_pc_contributing_genes The number of top genes contributing to
#'  each PC to include in the analysis. Defaults to 50.
#' @param experimental_design (Optional) A character string describing the 
#'  experimental design. This information helps the LLM contextualize the 
#'  analysis. If not provided, the LLM will analyze the data without specific
#'  experimental context.
#' @param temperature Controls the creativity of the LLM's response. Lower values 
#'  produce more deterministic results. Defaults to 0.
#' @param max_output_tokens The maximum number of tokens the LLM can generate. Defaults to 10048.
#' @param model_G Character string specifying the Gemini model to use for
#'  analysis. Default is "gemini-1.5-flash-latest".
#' @param model_O Character string specifying the Ollama model to use for
#'  analysis. Default is "llama3".
#' @param model_C Character string specifying the OpenAI model to use for
#'  analysis. Default is "gpt-4o-mini".
#' @param api_key_file The path to the file containing the API key for the LLM. Defaults to "api_keys.txt".
#'
#' @return A character string containing the LLM's analysis of the top PCs,
#'  including interpretations of the biological processes driving the variations,
#'  lists of top contributing genes and their associated pathways, and a summary
#'  paragraph.
#'
#' @usage SCassist_analyze_pcs(llm_server="google",
#'                    seurat_object_name, num_pcs = 5, 
#'                    top_n_pc_contributing_genes = 50, 
#'                    experimental_design = "",
#'                    temperature = 0,
#'                    max_output_tokens = 10048,
#'                    model_G = "gemini-1.5-flash-latest",
#'                    model_O = "llama3",
#'                    model_C = "gpt-4o-mini",
#'                    api_key_file = "api_keys.txt",
#'                    model_params = list(seed = 42, temperature = 0, num_gpu = 0))
#'
#' @import httr
#' @import jsonlite
#' @importFrom httr POST content content_type_json
#' @importFrom stats quantile
#' @importFrom jsonlite fromJSON
#' @keywords analyze, pca, principal components, large language model, llm
#' @export
SCassist_analyze_pcs <- function(llm_server="google",
                                     seurat_object_name,
                                     num_pcs = 5, 
                                     top_n_pc_contributing_genes = 50, 
                                     experimental_design = "",
                                     temperature = 0,
                                     max_output_tokens = 10048,
                                     model_G = "gemini-1.5-flash-latest",
                                     model_O = "llama3",
                                     model_C = "gpt-4o-mini",
                                     api_key_file = "api_keys.txt",
                                     model_params = list(seed = 42, temperature = 0, num_gpu = 0)
) {
  
  if (llm_server == "google") {
    return(SCassist_analyze_pcs_G(seurat_object_name,
                                      num_pcs = num_pcs, 
                                      top_n_pc_contributing_genes = top_n_pc_contributing_genes, 
                                      experimental_design = experimental_design,
                                      temperature = temperature,
                                      max_output_tokens = max_output_tokens,
                                      model = model_G,
                                      api_key_file = api_key_file ))
  } else if (llm_server == "ollama") {
    return(SCassist_analyze_pcs_L(seurat_object_name, 
                                      num_pcs = num_pcs, 
                                      top_n_pc_contributing_genes = top_n_pc_contributing_genes, 
                                      experimental_design = experimental_design,
                                      model_params = model_params, 
                                      model = model_O ))
  } else if (llm_server == "openai") {
    return(SCassist_analyze_pcs_C(seurat_object_name, 
                                      model_params = model_params, 
                                      num_pcs = num_pcs, 
                                      top_n_pc_contributing_genes = top_n_pc_contributing_genes, 
                                      experimental_design = experimental_design,
                                      temperature = temperature,
                                      max_output_tokens = max_output_tokens,
                                      model = model_C, 
                                      api_key_file = api_key_file))
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama'.")
  }
}

SCassist_analyze_pcs_G <- function(seurat_object_name, num_pcs = 5, 
                                 top_n_pc_contributing_genes = 50, 
                                 experimental_design = "",
                                 temperature = 0,
                                 max_output_tokens = 10048,
                                 model = "gemini-1.5-flash-latest",
                                 api_key_file = "api_keys.txt") {
  # Set up the Gemini API request
  model_query <- paste0(model, ":generateContent")  # Construct the API endpoint using the specified Gemini model name
  
  # Read the API key from the specified file
  api_key <- readLines(api_key_file)
  
  # Get the Seurat object
  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in environment.", call. = FALSE)
    }
  )
  
  # Check if PCA is available
  if (is.null(seurat_object@reductions$pca)) {
    stop("Please run PCA on the Seurat object first (e.g., using `RunPCA`).")
  }
  
  # Get the top genes contributing to each PC
  top_genes <- lapply(1:num_pcs, function(pc) {
    # Get gene loadings for the current PC
    gene_loadings <- seurat_object@reductions$pca@feature.loadings[, pc]
    # Get the top genes with highest absolute loadings
    top_genes <- names(sort(abs(gene_loadings), decreasing = TRUE))[1:top_n_pc_contributing_genes]
    return(top_genes)
  })
  
  # Initialize lists to store LLM responses and combined summary
  pc_summaries <- list()
  overall_summary <- ""
  
  # Iterate through each PC and analyze the top genes using LLM
  for (i in 1:num_pcs) {
    genes <- top_genes[[i]]
    
    # Create a prompt for the LLM for the current PC
    prompt <- paste0("PC", i, ": ", paste(genes, collapse = ", "), "\n\n")
    
    # Add instructions to the prompt
    full_prompt <- paste0(experimental_design,
                          "We performed QC, normalization and PCA on this data using Seurat. Here is the list of top PC's and their genes from our analysis:\n\n",
                          prompt, 
                          "\nIdentify the top contributing genes for this PC. Based on the gene functions and biological pathways associated with these genes, suggest potential biological processes that might be driving the variations captured by this PC. Present your results in a short paragraph\n\n")
    
    response1 <- POST(
      url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
      query = list(key = api_key), # Access api_key from api_keys.R
      content_type_json(),
      encode = "json",
      body = list(
        contents = list(
          parts = list(
            list(text = full_prompt)
          )),
        generationConfig = list(
          temperature = temperature,
          maxOutputTokens = max_output_tokens,
          seed = 123456  # For reproducible runs
        )
      )
    )
    
    candidates <- content(response1)$candidates
    outputs <- unlist(lapply(candidates, function(candidate) candidate$content$parts))
    
    # Store the LLM summary for this PC
    pc_summaries[[i]] <- outputs[["text"]]
    
    # Add to the overall summary 
    overall_summary <- paste0(overall_summary, "\n**PC", i, " Summary:**\n", outputs[["text"]], "\n")
  }
  
  # Combine all PC summaries into a single prompt
  combined_prompt <- paste0("Please provide an overall summary based on the individual summaries of each PC:\n", overall_summary)
  
  overall_response <- POST(
    url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
    query = list(key = api_key), # Access api_key from api_keys.R
    content_type_json(),
    encode = "json",
    body = list(
      contents = list(
        parts = list(
          list(text = combined_prompt)
        )),
      generationConfig = list(
        temperature = temperature,
        maxOutputTokens = max_output_tokens,
        seed = 123456  # For reproducible runs
      )
    )
  )
  
  candidates2 <- content(overall_response)$candidates
  outputs2 <- unlist(lapply(candidates2, function(candidate2) candidate2$content$parts))
  
  # Print the overall summary
  cat(overall_summary)
  cat("\n")
  cat(outputs2[["text"]])
  
  # Return the overall summary (for potential use elsewhere)
  return(outputs2[["text"]])
}

SCassist_analyze_pcs_L <- function(seurat_object_name, num_pcs = 5, 
                                   top_n_pc_contributing_genes = 50, 
                                   experimental_design = "",
                                   model_params = list(seed = 42, temperature = 0),
                                   model = "llama3") {
  # Get the Seurat object
  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in environment.", call. = FALSE)
    }
  )
  
  # Check if PCA is available
  if (is.null(seurat_object@reductions$pca)) {
    stop("Please run PCA on the Seurat object first (e.g., using `RunPCA`).")
  }
  
  # Get the top genes contributing to each PC
  top_genes <- lapply(1:num_pcs, function(pc) {
    # Get gene loadings for the current PC
    gene_loadings <- seurat_object@reductions$pca@feature.loadings[, pc]
    # Get the top genes with highest absolute loadings
    top_genes <- names(sort(abs(gene_loadings), decreasing = TRUE))[1:top_n_pc_contributing_genes]
    return(top_genes)
  })
  
  # Initialize lists to store LLM responses and combined summary
  pc_summaries <- list()
  overall_summary <- ""
  
  # Iterate through each PC and analyze the top genes using LLM
  for (i in 1:num_pcs) {
    genes <- top_genes[[i]]
    
    # Create a prompt for the LLM for the current PC
    prompt <- paste0("PC", i, ": ", paste(genes, collapse = ", "), "\n\n")
    
    # Add instructions to the prompt
    full_prompt <- paste0(experimental_design,
                          "We performed QC, normalization and PCA on this data using Seurat. Here is the list of top PC's and their genes from our analysis:\n\n",
                          prompt, 
                          "\nIdentify the top contributing genes for this PC. Based on the gene functions and biological pathways associated with these genes, suggest potential biological processes that might be driving the variations captured by this PC. Present your results in this example format:\n\n",
                          "PC# top genes.\n",
                          "Bulleted list of associated pathways.\n",
                          "One 3-line paragraph summarizing the result.\n")
    
    # Query the LLM
    response <- tryCatch(
      {
        suppressMessages(
          rollama::query(
            full_prompt,
            screen = FALSE,
            model_params = model_params,
            model = model
          )
        )
      },
      error = function(e) {
        stop("Error: The LLM query encountered an error. ",
             "Please check your ollama server connection and or the model.", call. = FALSE)
      }
    )
    
    # Check if LLM response is valid
    if (is.null(response[[1]]$message$content)) {
      stop("Error: The LLM returned an invalid response. Please check the LLM model and parameters.")
    }
    
    # Store the LLM summary for this PC
    pc_summaries[[i]] <- response[[1]]$message$content
    
    # Add to the overall summary 
    overall_summary <- paste0(overall_summary, "\n**PC", i, " Summary:**\n", response[[1]]$message$content, "\n")
  }
  
  # Combine all PC summaries into a single prompt
  combined_prompt <- paste0("Please provide an overall summary based on the individual summaries of each PC:\n", overall_summary)
  
  # Query the LLM for the overall summary
  overall_response <- tryCatch(
    {
      suppressMessages(
        rollama::query(
          combined_prompt,
          screen = FALSE,
          model_params = model_params,
          model = model
        )
      )
    },
    error = function(e) {
      stop("Error: The LLM query encountered an error. ",
           "Please check your ollama server connection and or the model.", call. = FALSE)
    }
  )
  
  # Check if LLM response is valid
  if (is.null(overall_response[[1]]$message$content)) {
    stop("Error: The LLM returned an invalid response. Please check the LLM model and parameters.")
  }
  
  # Print the overall summary
  cat(overall_summary)
  cat("\n")
  cat(overall_response[[1]]$message$content)
  
  # Return the overall summary (for potential use elsewhere)
  return(overall_response[[1]]$message$content)
}


SCassist_analyze_pcs_C <- function(seurat_object_name, num_pcs = 5, 
                                   top_n_pc_contributing_genes = 50, 
                                   experimental_design = "",
                                   temperature = 0,
                                   max_output_tokens = 10048,
                                   model = "gpt-4o-mini",
                                   model_params = list(seed = 42, temperature = 0),
                                   api_key_file = api_key_file) {
  # 1. Read the API key from the specified file
  api_key <- readLines(api_key_file, encoding = "UTF-8")
  
  # 2. Retrieve the Seurat object
  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in environment.", call. = FALSE)
    }
  )
  
  # Check if PCA is available
  if (is.null(seurat_object@reductions$pca)) {
    stop("Please run PCA on the Seurat object first (e.g., using `RunPCA`).")
  }
  
  # Get the top genes contributing to each PC
  top_genes <- lapply(1:num_pcs, function(pc) {
    # Get gene loadings for the current PC
    gene_loadings <- seurat_object@reductions$pca@feature.loadings[, pc]
    # Get the top genes with highest absolute loadings
    top_genes <- names(sort(abs(gene_loadings), decreasing = TRUE))[1:top_n_pc_contributing_genes]
    return(top_genes)
  })
  
  # Initialize lists to store LLM responses and combined summary
  pc_summaries <- list()
  overall_summary <- ""
  
  # Iterate through each PC and analyze the top genes using LLM
  for (i in 1:num_pcs) {
    genes <- top_genes[[i]]
    
    # Create a prompt for the LLM for the current PC
    prompt <- paste0("PC", i, ": ", paste(genes, collapse = ", "), "\n\n")
    
    # Add instructions to the prompt
    full_prompt <- paste0(experimental_design,
                          "We performed QC, normalization and PCA on this data using Seurat. Here is the list of top PC's and their genes from our analysis:\n\n",
                          prompt, 
                          "\nIdentify the top contributing genes for this PC. Based on the gene functions and biological pathways associated with these genes, suggest potential biological processes that might be driving the variations captured by this PC. Present your results in a short paragraph\n\n")

    # 11. Send the request to the OpenAI API
    # Make the POST request to OpenAI
    response1 <- httr::POST(
      url = "https://api.openai.com/v1/chat/completions",
      add_headers("Authorization" = paste("Bearer", api_key)),
      content_type_json(), 
      encode = "json",
      body = list(
        model = model, # Specify the model
        messages = list(      # Messages array as required by Chat Completions API
          list(
            role = "user",
            content = full_prompt
          )
        )
      )
    )
    
    # 12. Extract the OpenAI model's response
    # Check for HTTP errors
    if (http_error(response1)) {
      stop(paste("OpenAI API request failed with status", http_status(response1)$message, "\nContent:", content(response1, "text")))
    }
    
    # Parse the JSON response
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Package 'jsonlite' is required but not installed.")
    }
    
    response_content <- content(response1, "text")
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
    
    # Store the LLM summary for this PC
    pc_summaries[[i]] <- generated_text
    
    # Add to the overall summary 
    overall_summary <- paste0(overall_summary, "\n**PC", i, " Summary:**\n", generated_text, "\n")
  }
  
  # Combine all PC summaries into a single prompt
  combined_prompt <- paste0("Please provide an overall summary based on the individual summaries of each PC:\n", overall_summary)
  
  # Make the POST request to OpenAI
  overall_response <- httr::POST(
    url = "https://api.openai.com/v1/chat/completions",
    add_headers("Authorization" = paste("Bearer", api_key)),
    content_type_json(), 
    encode = "json",
    body = list(
      model = model, # Specify the model
      messages = list(      # Messages array as required by Chat Completions API
        list(
          role = "user",
          content = combined_prompt
        )
      )
    )
  )
  
  # 12. Extract the OpenAI model's response
  # Check for HTTP errors
  if (http_error(overall_response)) {
    stop(paste("OpenAI API request failed with status", http_status(overall_response)$message, "\nContent:", content(overall_response, "text")))
  }
  
  # Parse the JSON response
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required but not installed.")
  }
  
  response_content <- content(overall_response, "text")
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
  
  outputs2 <- generated_text
  
  # Print the overall summary
  cat(overall_summary)
  cat("\n")
  cat(outputs2)
  
  # Return the overall summary (for potential use elsewhere)
  return(outputs2)
}