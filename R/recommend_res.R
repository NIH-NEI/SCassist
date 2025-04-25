#' @title Recommend Resolution for Seurat's `FindClusters` Function
#'
#' @description This function analyzes key characteristics of a single-cell RNA-seq dataset
#' using a large language model (LLM) to recommend a suitable resolution range
#' for Seurat's `FindClusters` function. The LLM considers factors like the
#' number of cells, the number of highly variable genes (HVGs), the mean
#' expression variability across genes, and the median neighbor distance in the
#' k-nearest neighbor graph to provide a reasoned suggestion.
#'
#' @author Vijay Nagarajan, PhD, NEI/NIH
#' 
#' @details This function was written with assistance from Google's Gemini and Meta's Llama3.
#' 
#' @param llm_server The LLM server to use. Options are "google" or "ollama" or "openai". Default is "google".
#' @param model_params A list of parameters to be passed to the `ollama::query` function.
#'   This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.
#' @param seurat_object_name The name of the Seurat object containing the
#'  single-cell RNA-seq data. The object should be accessible in the current
#'  environment and should have `FindNeighbors` already run with 
#'  `return.neighbor = TRUE`.
#' @param use_SCT  Logical value indicating whether to use SCT neighbor distances
#'  (TRUE) or RNA neighbor distances (FALSE). Defaults to **TRUE** (SCT). If you
#'  ran SCTransform
#' @param temperature Numeric value between 0 and 1, controlling the
#'   "creativity" of the Gemini model's response. Higher values lead to more
#'   diverse and potentially unexpected outputs. Default is 0.
#' @param max_output_tokens Integer specifying the maximum number of tokens
#'   the Gemini model can generate in its response. Default is 10048.
#' @param model_G Character string specifying the Gemini model to use for
#'  analysis. Default is "gemini-1.5-flash-latest".
#' @param model_O Character string specifying the Ollama model to use for
#'  analysis. Default is "llama3".
#' @param model_C Character string specifying the OpenAI model to use for
#'  analysis. Default is "gpt-4o-mini".
#' @param api_key_file Character string specifying the path to a file containing
#'   the API key for accessing the Gemini model.
#'
#' @return A character string containing the LLM's recommendation for a 
#'  resolution range, along with a justification for the choice. The output
#'  will be structured as:
#'  - "Recommended Resolution:" followed by a range of values.
#'  - "Reasoning:" followed by a paragraph explaining the rationale behind
#'  the recommendation.
#'
#' @usage SCassist_recommend_res(llm_server="google",
#'        seurat_object_name, use_SCT = TRUE,
#'        temperature = 0,
#'        max_output_tokens = 10048,
#'        model_G = "gemini-1.5-flash-latest",
#'        model_O = "llama3",
#'        model_C = "gpt-4o-mini",
#'        api_key_file = "api_keys.txt",
#'        model_params = list(seed = 42, temperature = 0, num_gpu = 0))
#' 
#' @examples
#' \dontrun{
#' # Assuming you have a Seurat object named 'seurat_obj' with FindNeighbors run
#' SCassist_recommend_res(seurat_object_name = "seurat_obj",
#'                        use_SCT = TRUE,
#'                        api_key_file = "my_api_key.txt")
#' }
#' @import rollama
#' @import jsonlite
#' @importFrom httr POST content content_type_json
#' @importFrom jsonlite fromJSON
#' @importFrom stats median var
#' 
#' @keywords single-cell RNA-seq, Seurat, FindClusters, resolution, LLM, recommendation
#' @export

SCassist_recommend_res <- function(llm_server="google",
                                     seurat_object_name,
                                     use_SCT = TRUE,
                                     temperature = 0,
                                     max_output_tokens = 10048,
                                     model_G = "gemini-1.5-flash-latest",
                                     model_O = "llama3",
                                     model_C = "gpt-4o-mini",
                                     api_key_file = "api_keys.txt",
                                     model_params = list(seed = 42, temperature = 0, num_gpu = 0)
) {
  
  if (llm_server == "google") {
    return(SCassist_recommend_res_G(seurat_object_name,
                                      use_SCT = use_SCT,
                                      temperature = temperature,
                                      max_output_tokens = max_output_tokens,
                                      model = model_G,
                                      api_key_file = api_key_file ))
  } else if (llm_server == "ollama") {
    return(SCassist_recommend_res_L(seurat_object_name, 
                                      use_SCT = use_SCT,
                                      model_params = model_params, 
                                      model = model_O ))
  } else if (llm_server == "openai") {
    return(SCassist_recommend_res_C(seurat_object_name, 
                                    use_SCT = use_SCT,
                                    temperature = temperature,
                                    max_output_tokens = max_output_tokens,
                                    model_params = model_params, 
                                    model = model_C,
                                    api_key_file = api_key_file))
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama'.")
  }
}

SCassist_recommend_res_G <- function(seurat_object_name, use_SCT = TRUE,
                                   temperature = 0,
                                   max_output_tokens = 10048,
                                   model = "gemini-1.5-flash-latest",
                                   api_key_file = "api_keys.txt") { 
  
  # 1. Set up the Gemini API request
  model_query <- paste0(model, ":generateContent")
  
  # Read the API key from the specified file
  api_key <- readLines(api_key_file)
  
  # 2. Get the Seurat object
  seurat_object <- tryCatch(
    get(seurat_object_name),
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in the environment. ", 
           "Make sure the object exists and is accessible.", call. = FALSE)
    }
  )
  
  # 3. Check if FindNeighbors is run and has neighbor distances
  if (use_SCT) {
    if (is.null(seurat_object@neighbors$SCT.nn@nn.dist)) {
      stop("Error: FindNeighbors must be run with return.neighbor = TRUE option for SCT. ",
           "Please rerun FindNeighbors with this option after running SCTransform.")
    }
  } else {
    if (is.null(seurat_object@neighbors$RNA.nn@nn.dist)) {
      stop("Error: FindNeighbors must be run with return.neighbor = TRUE option for RNA. ",
           "Please rerun FindNeighbors with this option.")
    }
  }
  
  # 4. Get basic data about your Seurat object
  n_cells <- ncol(seurat_object)
  n_genes <- nrow(seurat_object)
  
  # 5. Calculate additional heterogeneity measures (examples):
  # - Mean expression variability:
  mean_expression_variability <- mean(apply(seurat_object@assays$RNA$counts, 1, var))
  
  # - Median neighbor distance:
  if (use_SCT) {
    median_neighbor_distance <- median(seurat_object@neighbors$SCT.nn@nn.dist)
  } else {
    median_neighbor_distance <- median(seurat_object@neighbors$RNA.nn@nn.dist)
  }
  
  # 6. Prepare data for LLM input
  data_for_llm <- list(
    n_cells = n_cells,
    n_genes = n_genes,
    #n_hvg = n_hvg, # Now based on the variance threshold
    mean_expression_variability = mean_expression_variability,
    median_neighbor_distance = median_neighbor_distance
  )
  
  
  # 7. LLM Prompt (using the data_for_llm variables)
  prompt <- paste0(
    "I am analyzing single-cell RNA sequencing data. My dataset consists of ",
    data_for_llm$n_cells, " cells and ", data_for_llm$n_genes, " genes. \n",
    "The mean expression variability across genes is ", data_for_llm$mean_expression_variability, 
    "and the median neighbor distance in the k-nearest neighbor graph is ", 
    data_for_llm$median_neighbor_distance, ". \n\n",
    "What resolution range would be most suitable for identifying distinct and subtle populations of cells in my data? ",
    "Provide the output as, Recommended Resolution, EXAMPLE: seq(starting resolution number,ending resolution number,increment number), : and a short reasoning paragraph under Reasoning. Start your response saying that, based on the data characterisitcs i recommend; "
  )
  
  # 8. Send the request to the Gemini API
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
  
  # 9. Extract the Gemini model's response
  candidates <- content(response1)$candidates
  outputs <- unlist(lapply(candidates, function(candidate) candidate$content$parts))
  
  # 10. Print and return the LLM's response
  cat(outputs[["text"]])
  return(outputs[["text"]])
}

SCassist_recommend_res_L <- function(seurat_object_name, use_SCT = TRUE,
                                     model_params = list(seed = 42, temperature = 0),
                                     model = "llama3") { 
  # Get the Seurat object
  seurat_object <- tryCatch(
    get(seurat_object_name),
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in the environment. ", 
           "Make sure the object exists and is accessible.", call. = FALSE)
    }
  )
  
  # Check if FindNeighbors is run and has neighbor distances
  if (use_SCT) {
    if (is.null(seurat_object@neighbors$SCT.nn@nn.dist)) {
      stop("Error: FindNeighbors must be run with return.neighbor = TRUE option for SCT. ",
           "Please rerun FindNeighbors with this option after running SCTransform.")
    }
  } else {
    if (is.null(seurat_object@neighbors$RNA.nn@nn.dist)) {
      stop("Error: FindNeighbors must be run with return.neighbor = TRUE option for RNA. ",
           "Please rerun FindNeighbors with this option.")
    }
  }
  
  # Get basic data about your Seurat object
  n_cells <- ncol(seurat_object)
  n_genes <- nrow(seurat_object)
  
  # Calculate additional heterogeneity measures (examples):
  # - Mean expression variability:
  mean_expression_variability <- mean(apply(seurat_object@assays$RNA$counts, 1, var))
  
  # - Median neighbor distance:
  if (use_SCT) {
    median_neighbor_distance <- median(seurat_object@neighbors$SCT.nn@nn.dist)
  } else {
    median_neighbor_distance <- median(seurat_object@neighbors$RNA.nn@nn.dist)
  }
  
  # Prepare data for LLM input
  data_for_llm <- list(
    n_cells = n_cells,
    n_genes = n_genes,
    #n_hvg = n_hvg, # Now based on the variance threshold
    mean_expression_variability = mean_expression_variability,
    median_neighbor_distance = median_neighbor_distance
  )
  
  
  # LLM Prompt (using the data_for_llm variables)
  prompt <- paste0(
    "I am analyzing single-cell RNA sequencing data. My dataset consists of ",
    data_for_llm$n_cells, " cells and ", data_for_llm$n_genes, " genes. \n",
    "The mean expression variability across genes is ", data_for_llm$mean_expression_variability, 
    "and the median neighbor distance in the k-nearest neighbor graph is ", 
    data_for_llm$median_neighbor_distance, ". \n\n",
    "What resolution range would be most suitable for identifying distinct and subtle populations of cells in my data? ",
    "Provide the output as, Recommended Resolution, EXAMPLE: seq(starting resolution number,ending resolution number,increment number), : and a short reasoning paragraph under Reasoning. Start your response saying that, based on the data characterisitcs i recommend; "
  )
  
  # Query the LLM (using ollama) with error handling
  response <- tryCatch(
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
           "Please check your ollama server connection.", call. = FALSE)
    }
  )
  
  # Check if LLM response is valid
  if (is.null(response[[1]]$message$content)) {
    stop("Error: The LLM returned an invalid response. Please check the LLM model and parameters.")
  }
  
  # Print and return the LLM's response
  cat(response[[1]]$message$content) 
  return(response[[1]]$message$content)
}


SCassist_recommend_res_C <- function(seurat_object_name, use_SCT = TRUE,
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
  
  # 3. Check if FindNeighbors is run and has neighbor distances
  if (use_SCT) {
    if (is.null(seurat_object@neighbors$SCT.nn@nn.dist)) {
      stop("Error: FindNeighbors must be run with return.neighbor = TRUE option for SCT. ",
           "Please rerun FindNeighbors with this option after running SCTransform.")
    }
  } else {
    if (is.null(seurat_object@neighbors$RNA.nn@nn.dist)) {
      stop("Error: FindNeighbors must be run with return.neighbor = TRUE option for RNA. ",
           "Please rerun FindNeighbors with this option.")
    }
  }
  
  # 4. Get basic data about your Seurat object
  n_cells <- ncol(seurat_object)
  n_genes <- nrow(seurat_object)
  
  # 5. Calculate additional heterogeneity measures (examples):
  # - Mean expression variability:
  mean_expression_variability <- mean(apply(seurat_object@assays$RNA$counts, 1, var))
  
  # - Median neighbor distance:
  if (use_SCT) {
    median_neighbor_distance <- median(seurat_object@neighbors$SCT.nn@nn.dist)
  } else {
    median_neighbor_distance <- median(seurat_object@neighbors$RNA.nn@nn.dist)
  }
  
  # 6. Prepare data for LLM input
  data_for_llm <- list(
    n_cells = n_cells,
    n_genes = n_genes,
    #n_hvg = n_hvg, # Now based on the variance threshold
    mean_expression_variability = mean_expression_variability,
    median_neighbor_distance = median_neighbor_distance
  )
  
  
  # 7. LLM Prompt (using the data_for_llm variables)
  prompt <- paste0(
    "I am analyzing single-cell RNA sequencing data. My dataset consists of ",
    data_for_llm$n_cells, " cells and ", data_for_llm$n_genes, " genes. \n",
    "The mean expression variability across genes is ", data_for_llm$mean_expression_variability, 
    "and the median neighbor distance in the k-nearest neighbor graph is ", 
    data_for_llm$median_neighbor_distance, ". \n\n",
    "What resolution range would be most suitable for identifying distinct and subtle populations of cells in my data? ",
    "Provide the output as, Recommended Resolution, EXAMPLE: seq(starting resolution number,ending resolution number,increment number), : and a short reasoning paragraph under Reasoning. Start your response saying that, based on the data characterisitcs i recommend; "
  )
  
  # 8. Send the request to the Gemini API
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