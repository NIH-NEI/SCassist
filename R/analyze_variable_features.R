#' @title Analyze Variable Features with a Large Language Model (LLM)
#' @description This function analyzes a list of variable features from a Seurat object
#' using a large language model (LLM) to identify enriched gene ontologies or
#' pathways. The LLM considers the provided gene list and experimental design
#' information to provide an analysis of the functional significance of these
#' variable features.
#' @author Vijay Nagarajan, PhD, NEI/NIH
#' @details This function was written with assistance from Google's Gemini and Meta's Llama3.
#' @seealso  \code{\link{SCassist_analyze_quality}}, 
#'           \code{\link{SCassist_recommend_normalization}}, 
#'           \code{\link{SCassist_analyze_variable_features}}
#' @usage 
#' SCassist_analyze_variable_features(llm_server = "google",
#'                  seurat_object_name, 
#'                  top_n_variable_features = 30, 
#'                  experimental_design = NULL, 
#'                  temperature = 0,
#'                  max_output_tokens = 10048,
#'                  model_G = "gemini-1.5-flash-latest",
#'                  model_O = "llama3",
#'                  model_C = "gpt-4o-mini",
#'                  api_key_file = "api_keys.txt",
#'                  model_params = list(seed = 42,temperature = 0, num_gpu = 0))
#' @param llm_server The LLM server to use. Options are "google" or "ollama" or "openai". Default is "google".
#' @param seurat_object_name The name of the Seurat object containing the
#'  single-cell RNA-seq data. The object should be accessible in the current
#'  environment and should have either `scTransform` or `FindVariableFeatures`
#'  already run.
#' @param top_n_variable_features The number of top variable features to analyze. 
#'  Defaults to 30 if not provided.
#' @param experimental_design A character string describing the experimental
#'  design. This information helps the LLM interpret the significance of
#'  the identified gene ontologies or pathways.  If not provided, the LLM will
#'  analyze the variable features without specific design context.
#' @param temperature A number between 0 and 1 controlling the creativity of the
#'  LLM's response. Lower values produce more deterministic results. Defaults
#'  to 0.
#' @param max_output_tokens The maximum number of tokens the LLM can generate 
#'  in its response. Defaults to 10048.
#' @param model_G Character string specifying the Gemini model to use for
#'  analysis. Default is "gemini-1.5-flash-latest".
#' @param model_O Character string specifying the Ollama model to use for
#'  analysis. Default is "llama3".
#' @param model_C Character string specifying the OpenAI model to use for
#'  analysis. Default is "gpt-4o-mini".
#' @param model_params A list of parameters to be passed to the `ollama::query` function.
#'   This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.
#' @param api_key_file The path to a text file containing the API key for 
#'  accessing the LLM.
#'
#' @return A character string containing the LLM's analysis of the variable
#'  features, including a table with enriched functional categories and scores,
#'  and an explanation of the relevance to the experimental design (if provided).
#'
#' @importFrom Seurat VariableFeatures
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the Seurat object
#' seurat_object <- readRDS("path/to/seurat_object.rds")
#' 
#' # Analyze variable features with Gemini
#' analysis_results <- SCassist_analyze_variable_features(
#'   seurat_object_name = "seurat_object",
#'   top_n_variable_features = 50,
#'   experimental_design = "Treatment with drug X compared to control",
#'   model = "gemini-1.5-flash-latest",
#'   api_key_file = "api_keys.txt"
#' )
#' 
#' # Print the analysis results
#' cat(analysis_results)
#' }
#' @import jsonlite
#' @importFrom httr POST content content_type_json
#' @export

SCassist_analyze_variable_features <- function(llm_server="google",
                                     seurat_object_name,
                                     top_n_variable_features = 30, 
                                     experimental_design = NULL, 
                                     temperature = 0,
                                     max_output_tokens = 10048,
                                     model_G = "gemini-1.5-flash-latest",
                                     model_O = "llama3",
                                     model_C = "gpt-4o-mini",
                                     api_key_file = "api_keys.txt",
                                     model_params = list(seed = 42, temperature = 0, num_gpu = 0)
) {
  
  if (llm_server == "google") {
    return(SCassist_analyze_variable_features_G(seurat_object_name,
                                      top_n_variable_features = top_n_variable_features, 
                                      experimental_design = experimental_design, 
                                      temperature = temperature,
                                      max_output_tokens = max_output_tokens,
                                      model = model_G,
                                      api_key_file = api_key_file ))
  } else if (llm_server == "ollama") {
    return(SCassist_analyze_variable_features_L(seurat_object_name, 
                                      top_n_variable_features = top_n_variable_features, 
                                      experimental_design = experimental_design,
                                      model_params = model_params, 
                                      model = model_O ))
  } else if (llm_server == "openai") {
    return(SCassist_analyze_variable_features_C(seurat_object_name, 
                                      top_n_variable_features = top_n_variable_features, 
                                      experimental_design = experimental_design,
                                      model_params = model_params, 
                                      model = model_C, 
                                      api_key_file = api_key_file))
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama'.")
  }
}

SCassist_analyze_variable_features_G <- function(seurat_object_name, 
                                               top_n_variable_features = 30, 
                                               experimental_design = NULL, 
                                               temperature = 0,
                                               max_output_tokens = 10048,
                                               model = "gemini-1.5-flash-latest",
                                               api_key_file = "api_keys.txt" ) {
  
  # 1. Set up the Gemini API request
  model_query <- paste0(model, ":generateContent")
  
  # 2. Read the API key from the specified file
  api_key <- readLines(api_key_file)
  
  # 3. Get the Seurat object
  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in environment.", call. = FALSE)
    }
  )
  
  # 4. Check if variable features are available
  if (is.null(VariableFeatures(seurat_object))) {
    stop("Please run `scTransform` or `FindVariableFeatures` on the Seurat object first.")
  }
  
  # 5. Prepare the Variable Feature Data
  variable_features <- head(VariableFeatures(object = seurat_object), top_n_variable_features)
  gene_list <- paste(variable_features, collapse = ", ")
  
  # 6. Craft a Prompt for the LLM
  prompt <- paste0("Analyze the following list of genes: ", gene_list,
                   "These genes were identified as variable features. Identify the most enriched gene ontologies or pathways among these genes. Dont provide any pvalues or enrichment scores, just summarize based on these genes known functions. Dont provide any Further Analysis suggestions")
  
  # 7. Add experimental design if provided
  if (!is.null(experimental_design)) {
    prompt <- paste0(prompt, "  Explain the relevance of these categories to the experimental design: ", experimental_design, ".")
  }
  
  # 8. Complete the prompt
  #prompt <- paste0(prompt, "  Provide your answer as a table with two columns: 'Functional Category' and 'Enrichment Score'.")
  
  # 9. Send the request to the Gemini API
  response1 <- POST(
    url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
    query = list(key = api_key), # Access api_key from api_keys.R
    content_type_json(),
    encode = "json",
    body = list(
      contents = list(
        parts = list(
          list(text = prompt)
        )),
      generationConfig = list(
        temperature = temperature,
        maxOutputTokens = max_output_tokens,
        seed = 123456  # For reproducible runs
      )
    )
  )
  
  # 10. Extract the Gemini model's response
  candidates <- content(response1)$candidates
  outputs <- unlist(lapply(candidates, function(candidate) candidate$content$parts))
  
  # 11. Print and return the LLM's response
  cat(outputs[["text"]])
  return(outputs[["text"]])
}

SCassist_analyze_variable_features_L <- function(seurat_object_name, top_n_variable_features = 30, 
                                                 experimental_design = NULL, model_params = list(seed = 42, temperature = 0), 
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
  
  # Check if variable features are available
  if (is.null(VariableFeatures(seurat_object))) {
    stop("Please run `scTransform` or `FindVariableFeatures` on the Seurat object first.")
  }
  
  # Prepare the Variable Feature Data
  variable_features <- head(VariableFeatures(object = seurat_object), top_n_variable_features)
  gene_list <- paste(variable_features, collapse = ", ")
  
  # Craft a Prompt for the LLM
  prompt <- paste0("Analyze the following list of genes: ", gene_list,
                   "These genes were identified as variable features. Identify the most enriched gene ontologies or pathways among these genes. Dont provide any pvalues or enrichment scores, just summarize based on these genes known functions.")
  
  # Add experimental design if provided
  if (!is.null(experimental_design)) {
    prompt <- paste0(prompt, "  Explain the relevance of these categories to the experimental design: ", experimental_design, ".")
  }
  
  # Complete the prompt
  #prompt <- paste0(prompt, "  Provide your answer as a table with two columns: 'Functional Category' and 'Enrichment Score'.")
  
  # Query the LLM
  response_vf <- tryCatch(
    {
      suppressMessages(
        rollama::query(
          prompt,
          model = model,
          screen = FALSE,
          model_params = model_params
        )
      )
    },
    error = function(e) {
      stop("Error: The LLM query encountered an error. ",
           "Please check your ollama server connection and/or the model.", call. = FALSE)
    }
  )
  
  # Check if LLM response is valid
  if (is.null(response_vf[[1]]$message$content)) {
    stop("Error: The LLM returned an invalid response. Please check the LLM model and parameters.")
  }
  
  # Print and return the LLM's response
  cat(response_vf[[1]]$message$content) 
  return(response_vf[[1]]$message$content)
}


SCassist_analyze_variable_features_C <- function(seurat_object_name, 
                                                 top_n_variable_features = 30, 
                                                 experimental_design = NULL, 
                                                 temperature = 0,
                                                 max_output_tokens = 10048,
                                                 model = "gpt-4o-mini",
                                                 model_params = model_params,
                                                 api_key_file = "api_keys.txt" ) {
  
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
  
  # 4. Check if variable features are available
  if (is.null(VariableFeatures(seurat_object))) {
    stop("Please run `scTransform` or `FindVariableFeatures` on the Seurat object first.")
  }
  
  # 5. Prepare the Variable Feature Data
  variable_features <- head(VariableFeatures(object = seurat_object), top_n_variable_features)
  gene_list <- paste(variable_features, collapse = ", ")
  
  # 6. Craft a Prompt for the LLM
  prompt <- paste0("Analyze the following list of genes: ", gene_list,
                   "These genes were identified as variable features. Identify the most enriched gene ontologies or pathways among these genes. Dont provide any pvalues or enrichment scores, just summarize based on these genes known functions. Dont provide any Further Analysis suggestions")
  
  # 7. Add experimental design if provided
  if (!is.null(experimental_design)) {
    prompt <- paste0(prompt, "  Explain the relevance of these categories to the experimental design: ", experimental_design, ".")
  }
  
  # 8. Complete the prompt
  #prompt <- paste0(prompt, "  Provide your answer as a table with two columns: 'Functional Category' and 'Enrichment Score'.")
  
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
          content = prompt
        )
      )
    )
  )
  
  # 12. Extract the OpenAI model's response
  # Check for HTTP errors
  if (http_error(response)) {
    stop(paste("OpenAI API request failed with status", http_status(response)$message, "\nContent:", content(response, "text")))
  }
  
  # Parse the JSON response
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required but not installed.")
  }
  
  response_content <- content(response, "text")
  response_json <- jsonlite::fromJSON(response_content, flatten = TRUE)
  
  
  # Extract and print the generated text from the response
  # For chat completions, the text is in choices[[1]]$message$content
  if (!is.null(response_json$choices) && length(response_json$choices) > 0 && !is.null(response_json$choices$message.content)) {
    generated_text <- response_json$choices$message.content
    cat(generated_text)
  } else {
    cat("Error: Could not extract generated text from OpenAI response.\n")
    print(response_json) # Print the full response for debugging
  }
  
  # 13. Print and return the LLM's response
  #cat(generated_text)
  return(generated_text)
}