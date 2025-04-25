#' @title Analyze Top Markers, Predict Cell Types, and Optionally Annotate a Seurat Object
#' @description This function performs a comprehensive analysis of single-cell data, including:
#' 1. **Analyzing top markers for each cluster:** It identifies the most highly expressed genes (markers) in each cluster of cells.
#' 2. **Predicting potential cell types:** It uses a large language model (LLM) to analyze the top markers and predict the cell type that is most likely represented by each cluster.
#' 3. **Optionally annotating a Seurat object:** If you provide the name of a Seurat object, the function will add a new column to the object's metadata containing the predicted cell types for each cell.
#'
#' **Important Note:** The function assumes that the current identity of the Seurat
#' object (set with `Idents(seurat_object)`) is the same as the identity used
#' when running `FindAllMarkers`. This ensures correct cluster annotation.
#' @author Vijay Nagarajan, PhD, NEI/NIH
#' @details This function was written with assistance from Google's Gemini.
#' @param all_markers A data frame containing markers for each cluster, as output
#'  by the `FindAllMarkers` function from the Seurat package. This data frame should
#'  have columns for cluster, gene, avg_logFC, p_val, and pct.1.
#'
#' @param llm_server The LLM server to use. Options are "google" or "ollama" or "openai". Default is "google".
#' @param model_params A list of parameters to be passed to the `ollama::query` function.
#'   This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.
#' @param seurat_object_name Character string representing the name of the
#'   Seurat object in the current environment.
#' @param top_genes The number of top markers to consider for each cluster. This
#'  parameter controls the number of genes used as input for the LLM. A higher
#'  value will provide more information to the LLM, but may also increase the
#'  processing time.
#'
#' @param temperature A value between 0 and 1 that controls the creativity of the LLM's response.
#'  Lower values produce more deterministic results, while higher values allow for
#'  more diverse and potentially unexpected predictions.
#'
#' @param max_output_tokens The maximum number of tokens the LLM can generate in its response.
#'  This parameter limits the length of the LLM's output. A higher value allows for
#'  more detailed explanations, but may also increase the processing time.
#'
#' @param model_G Character string specifying the Gemini model to use for
#'  analysis. Default is "gemini-1.5-flash-latest".
#'
#' @param model_O Character string specifying the Ollama model to use for
#'  analysis. Default is "llama3".
#'
#' @param model_C Character string specifying the OpenAI model to use for
#'  analysis. Default is "gpt-4o-mini".
#'  
#' @param api_key_file The path to a file containing your Gemini API key. This key
#'  is required to access the Gemini API. You can obtain a key from Google Cloud
#'  Platform.
#'
#' @param seurat_object_name (optional) The name of the Seurat object to annotate.
#'  If provided, the function will annotate the object with predicted cell types.
#'  Default is NULL, which skips annotation.
#'
#' @return Returns a data frame containing the predicted cell types and reasoning
#'  for each cluster (`sca_annotation`). The data frame has three columns:
#'  - `cluster_num`: The cluster number.
#'  - `cell_type`: The predicted cell type.
#'  - `scassist_reasoning`: A brief explanation of the reasoning behind the prediction.
#'
#'  The function also updates the Seurat object in place if `seurat_object_name` is provided.
#'  The updated object will have a new column named `scassist_annotation_` followed by the
#'  name of the current identity column, containing the predicted cell types for each cell.
#'
#' @usage
#' SCassist_analyze_and_annotate(llm_server="google",
#'                              seurat_object_name = NULL,
#'                              all_markers,
#'                              top_genes = 30,
#'                              temperature = 0,
#'                              max_output_tokens = 10048,
#'                              model_G = "gemini-1.5-flash-latest",
#'                              model_O = "llama3",
#'                              model_C = "gpt-4o-mini",
#'                              api_key_file = "api_keys.txt",
#'                              model_params = list(seed = 42,
#'                                    temperature = 0, num_gpu = 0)
#'                              )
#'
#' @examples
#' \dontrun{
#' # Assuming you have a Seurat object named 'seurat_obj' and a data frame named 'markers'
#' # containing the results of FindAllMarkers
#' SCassist_analyze_and_annotate(all_markers = markers,
#'                              seurat_object_name = "seurat_obj",
#'                              api_key_file = "my_api_key.txt")
#' }
#' @keywords single-cell analysis, cell type prediction, Seurat, LLM, annotation
#' @import httr
#' @import jsonlite
#' @importFrom httr POST content content_type_json
#' @importFrom jsonlite fromJSON
#' @importFrom utils write.table
#' @importFrom plyr mapvalues
#' @importFrom Seurat Idents
#' @export

SCassist_analyze_and_annotate <- function(llm_server="google",
                                     seurat_object_name = NULL,
                                     all_markers, top_genes = 30,
                                     temperature = 0,
                                     max_output_tokens = 10048,
                                     model_G = "gemini-1.5-flash-latest",
                                     model_O = "llama3",
                                     model_C = "gpt-4o-mini",
                                     api_key_file = "api_keys.txt",
                                     model_params = list(seed = 42, temperature = 0, num_gpu = 0)
) {

  if (llm_server == "google") {
    return(SCassist_analyze_and_annotate_G(seurat_object_name = seurat_object_name,
                                      all_markers = all_markers,
                                      top_genes = top_genes,
                                      temperature = temperature,
                                      max_output_tokens = max_output_tokens,
                                      model = model_G,
                                      api_key_file = api_key_file ))
  } else if (llm_server == "ollama") {
    return(SCassist_analyze_and_annotate_L(seurat_object_name = seurat_object_name,
                                      all_markers = all_markers,
                                      top_genes = top_genes,
                                      model_params = model_params,
                                      model = model_O ))
  } else if (llm_server == "openai") {
    return(SCassist_analyze_and_annotate_C(seurat_object_name = seurat_object_name,
                                           all_markers = all_markers,
                                           top_genes = top_genes,
                                           temperature = temperature,
                                           max_output_tokens = max_output_tokens,
                                           model_params = model_params,
                                           model = model_C,
                                           api_key_file = api_key_file))
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama'.")
  }
}

SCassist_analyze_and_annotate_G <- function(all_markers, top_genes = 30,
                                          temperature = 0,
                                          max_output_tokens = 10048,
                                          model = "gemini-1.5-flash-latest",
                                          api_key_file = "api_key.txt",
                                          seurat_object_name = NULL) {

  # 1. Set up the Gemini API request
  model_query <- paste0(model, ":generateContent")  # Construct the API endpoint using the specified Gemini model name

  # Read the API key from the specified file
  api_key <- readLines(api_key_file)

  # Prompt template for the LLM
  prompt_template <- "The provided genes are the top markers of this single cell cluster. Analyze it and predict a potential cell type based on the markers. provide output in three columns. the first column should be the cluster number, second column should be the name of the potential cell type, third column should be a one paragraph reasoning. separate the columns using a colon. dont provide any other additional content generated by you, like: here is the analysis, etc. Here is an EXAMPLE output - '8:Megakaryocyte-precursor cells: The combination of markers LY6G6F, GP9, ITGA2B, and TMEM40 suggests a megakaryocytic origin, with involvement in platelet development and function'. Here is the input cluster number and corresponding markers for your analysis:"

  # 2. Prepare the input for the LLM
  top_markers_per_cluster <- split(all_markers, all_markers$cluster)  # Group markers by cluster
  top_markers_per_cluster <- lapply(top_markers_per_cluster, function(cluster_markers) {
    gene_names <- sub("\\..*", "", rownames(cluster_markers))  # Extract gene names
    gene_names[1:top_genes]  # Take the top 'top_genes' markers
  })

  # Printout message
  print("The output of Analyze and Annotate is saved as a tab delimited text if your current working directory. If you provided the name of the seurat object, the annotation is also added to that object in the SCassist_annotation column")

  # 3. Query the LLM for each cluster
  cluster_summaries <- list()  # Store LLM responses

  for (cluster_num in names(top_markers_per_cluster)) {
    cluster_markers <- top_markers_per_cluster[[cluster_num]]
    full_prompt <- paste0(prompt_template, "cluster ", cluster_num, ": ", paste(cluster_markers, collapse = ", "))  # Construct the prompt

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

    print(paste0("SCassistant is analyzing markers to predict cell types for cluster : ",cluster_num,collapse = ""))
    cluster_summaries[[cluster_num]] <- outputs[["text"]]  # Store the LLM response
  }

  # 4. Process the LLM responses
  sca_annotation <- data.frame(cluster_num = character(),
                               cell_type = character(),
                               scassist_reasoning = character(),
                               stringsAsFactors = FALSE)  # Create an empty data frame to store the annotation results

  for (i in seq_along(cluster_summaries)) {
    parts <- strsplit(cluster_summaries[[i]], ":")[[1]]  # Split the response into parts based on the colon delimiter

    cluster_num <- parts[1]
    cell_type <- parts[2]
    scassist_reasoning <- paste(parts[3], sep = "")  # Combine the reasoning parts

    cell_type <- trimws(cell_type)  # Remove leading and trailing spaces
    scassist_reasoning <- trimws(scassist_reasoning)

    sca_annotation <- rbind(sca_annotation, data.frame(cluster_num = cluster_num,
                                                       cell_type = cell_type,
                                                       scassist_reasoning = scassist_reasoning))
  }

  print(sca_annotation)  # Print the annotation results

  # 5. Save the annotation results
  suppressWarnings(
    write.table(sca_annotation, file = "G_sca_annotation_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  )

  # 6. Annotate the Seurat object (if requested)
  if (!is.null(seurat_object_name)) {
    seurat_object <- get(seurat_object_name, envir = parent.frame())  # Get the Seurat object from the parent environment

    current_identity <- Idents(seurat_object)  # Get the current identity of the Seurat object
    current_identity_column <- NULL  # Initialize the identity column name

    for (col_name in colnames(seurat_object@meta.data)) {
      if (all(as.character(current_identity) == as.character(seurat_object@meta.data[[col_name]]))) {
        current_identity_column <- col_name
        break
      }
    }

    if (is.null(current_identity_column)) {
      stop("Error: Unable to find the current identity column in the Seurat object. ",
           "Make sure the object has an identity set using 'Idents(seurat_object)'")
    }

    seurat_object[[paste0("scassist_annotation_", current_identity_column)]] <- tryCatch(
      {
        plyr::mapvalues(
          x = as.character(seurat_object@meta.data[[current_identity_column]]),  # Convert to character
          from = as.character(sca_annotation$cluster_num),  # Convert to character
          to = sca_annotation$cell_type
        )
      },
      error = function(e) {
        stop("Error: An error occurred while annotating clusters. ",
             "Please check if 'sca_annotation' is available and has the correct format.", call. = FALSE)
      }
    )
    # Return the modified Seurat object
    return(seurat_object)
    #assign(seurat_object_name, seurat_object, envir = parent.frame())  # Update the Seurat object in the parent environment
  }
  else {
    # Return the annotation results if no Seurat object is provided
    return(sca_annotation)
  }
  # 7. Return the annotation results
  #return(sca_annotation)
}

SCassist_analyze_and_annotate_L <- function(all_markers, top_genes = 30,
                                            model = "llama3",
                                            model_params = list(seed = 42, temperature = 0, num_gpu = 0),
                                            seurat_object_name = NULL) {

  # Prompt template
  prompt_template <- "The provided genes are the top markers of this single cell cluster. Analyze it and predict a potential cell type based on the markers. provide output in three columns. the first column should be the cluster number, second column should be the name of the potential cell type, third column should be a one paragraph reasoning. separate the columns using a colon. dont provide any other additional content generated by you, like: here is the analysis, etc. Here is an EXAMPLE output - '8:Megakaryocyte-precursor cells: The combination of markers LY6G6F, GP9, ITGA2B, and TMEM40 suggests a megakaryocytic origin, with involvement in platelet development and function'. Here is the input cluster number and corresponding markers for your analysis:"

  top_markers_per_cluster <- split(all_markers, all_markers$cluster)
  top_markers_per_cluster <- lapply(top_markers_per_cluster, function(cluster_markers) {
    gene_names <- sub("\\..*", "", rownames(cluster_markers))
    gene_names[1:top_genes]  # Extract the top 'top_genes' gene names
  })

  # Printout message
  print("The output of Analyze and Annotate is saved as a tab delimited text if your current working directory. If you provided the name of the seurat object, the annotation is also added to that object in the SCassist_annotation column")
  # Get top markers for each cluster

  # Store responses in a list
  cluster_summaries <- list()

  # Summarize top markers using Ollama for each cluster
  for (cluster_num in names(top_markers_per_cluster)) {
    cluster_markers <- top_markers_per_cluster[[cluster_num]]
    full_prompt <- paste0(prompt_template, "cluster ", cluster_num, ": ", paste(cluster_markers, collapse = ", "))

    response1 <- tryCatch(
      {
        suppressMessages(
          rollama::query(
            full_prompt,
            screen = FALSE,
            model = model,
            model_params = model_params
          )
        )
      },
      error = function(e) {
        stop("Error: The LLM query encountered an error. ",
             "Please check your ollama server connection.", call. = FALSE)
      }
    )

    print(paste0("SCassistant is analyzing markers to predict cell types for cluster : ",cluster_num,collapse = ""))
    cluster_summaries[[cluster_num]] <- response1[[1]]$message$content
  }

  # Create an empty data frame
  sca_annotation <- data.frame(cluster_num = character(),
                               cell_type = character(),
                               scassist_reasoning = character(),
                               stringsAsFactors = FALSE)

  # Loop through clusters and split the strings
  for (i in seq_along(cluster_summaries)) {
    parts <- strsplit(cluster_summaries[[i]], ":")[[1]]  # Split the string

    # Extract and combine the data
    cluster_num <- parts[1]
    cell_type <- parts[2]
    # Only take the first 3 parts, not 4
    scassist_reasoning <- paste(parts[3], sep = "") # Combine without space

    # Trim leading/trailing spaces
    cell_type <- trimws(cell_type)
    scassist_reasoning <- trimws(scassist_reasoning)

    # Add a new row to the data frame
    sca_annotation <- rbind(sca_annotation, data.frame(cluster_num = cluster_num,
                                                       cell_type = cell_type,
                                                       scassist_reasoning = scassist_reasoning))

    # Print the output to the screen (in addition to the file)
  }

  print(sca_annotation)

  # Write the sca_annotation data frame to a tab-separated file
  suppressWarnings(
    write.table(sca_annotation, file = "L_sca_annotation_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  )
  # Annotate the Seurat object if requested
  if (!is.null(seurat_object_name)) {
    # Get the Seurat object by name (using get)
    seurat_object <- get(seurat_object_name, envir = parent.frame())

    # Get the current identity column name
    current_identity <- Idents(seurat_object)
    current_identity_column <- NULL  # Initialize the column name

    # Loop through metadata columns and find a match
    for (col_name in colnames(seurat_object@meta.data)) {
      # Convert both to character vectors for comparison
      if (all(as.character(current_identity) == as.character(seurat_object@meta.data[[col_name]]))) {
        current_identity_column <- col_name
        break
      }
    }

    # Check if the identity column was found
    if (is.null(current_identity_column)) {
      stop("Error: Unable to find the current identity column in the Seurat object. ",
           "Make sure the object has an identity set using 'Idents(seurat_object)'")
    }

    seurat_object[[paste0("scassist_annotation_", current_identity_column)]] <- tryCatch(
      {
        plyr::mapvalues(
          x = as.character(seurat_object@meta.data[[current_identity_column]]),  # Convert to character
          from = as.character(sca_annotation$cluster_num),  # Convert to character
          to = sca_annotation$cell_type
        )
      },
      error = function(e) {
        stop("Error: An error occurred while annotating clusters. ",
             "Please check if 'sca_annotation' is available and has the correct format.", call. = FALSE)
      }
    )
    # Return the modified Seurat object
    return(seurat_object)
  } else {
    # Return the annotation results if no Seurat object is provided
    return(sca_annotation)
  }
}

SCassist_analyze_and_annotate_C <- function(all_markers, top_genes = 30,
                                            temperature = 0,
                                            max_output_tokens = 10048,
                                            model = "gpt-4o-mini",
                                            model_params = list(seed = 42, temperature = 0),
                                            api_key_file = api_key_file,
                                            seurat_object_name = NULL) {
  
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
  
  # Prompt template for the LLM
  prompt_template <- "The provided genes are the top markers of this single cell cluster. Analyze it and predict a potential cell type based on the markers. provide output in three columns. the first column should be the cluster number, second column should be the name of the potential cell type, third column should be a one paragraph reasoning. separate the columns using a colon. dont provide any other additional content generated by you, like: here is the analysis, etc. Here is an EXAMPLE output - '8:Megakaryocyte-precursor cells: The combination of markers LY6G6F, GP9, ITGA2B, and TMEM40 suggests a megakaryocytic origin, with involvement in platelet development and function'. Here is the input cluster number and corresponding markers for your analysis:"
  
  # 2. Prepare the input for the LLM
  top_markers_per_cluster <- split(all_markers, all_markers$cluster)  # Group markers by cluster
  top_markers_per_cluster <- lapply(top_markers_per_cluster, function(cluster_markers) {
    gene_names <- sub("\\..*", "", rownames(cluster_markers))  # Extract gene names
    gene_names[1:top_genes]  # Take the top 'top_genes' markers
  })
  
  # Printout message
  print("The output of Analyze and Annotate is saved as a tab delimited text if your current working directory. If you provided the name of the seurat object, the annotation is also added to that object in the SCassist_annotation column")
  
  # 3. Query the LLM for each cluster
  cluster_summaries <- list()  # Store LLM responses
  
  for (cluster_num in names(top_markers_per_cluster)) {
    cluster_markers <- top_markers_per_cluster[[cluster_num]]
    full_prompt <- paste0(prompt_template, "cluster ", cluster_num, ": ", paste(cluster_markers, collapse = ", "))  # Construct the prompt
    
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
            content = full_prompt
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
      #cat(generated_text)
    } else {
      cat("Error: Could not extract generated text from OpenAI response.\n")
      print(response_json) # Print the full response for debugging
    }
    
    print(paste0("SCassistant is analyzing markers to predict cell types for cluster : ",cluster_num,collapse = ""))
    cluster_summaries[[cluster_num]] <- generated_text  # Store the LLM response
  }
  
  # 4. Process the LLM responses
  sca_annotation <- data.frame(cluster_num = character(),
                               cell_type = character(),
                               scassist_reasoning = character(),
                               stringsAsFactors = FALSE)  # Create an empty data frame to store the annotation results
  
  for (i in seq_along(cluster_summaries)) {
    parts <- strsplit(cluster_summaries[[i]], ":")[[1]]  # Split the response into parts based on the colon delimiter
    
    cluster_num <- parts[1]
    cell_type <- parts[2]
    scassist_reasoning <- paste(parts[3], sep = "")  # Combine the reasoning parts
    
    cell_type <- trimws(cell_type)  # Remove leading and trailing spaces
    scassist_reasoning <- trimws(scassist_reasoning)
    
    sca_annotation <- rbind(sca_annotation, data.frame(cluster_num = cluster_num,
                                                       cell_type = cell_type,
                                                       scassist_reasoning = scassist_reasoning))
  }
  
  print(sca_annotation)  # Print the annotation results
  
  # 5. Save the annotation results
  suppressWarnings(
    write.table(sca_annotation, file = "C_sca_annotation_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  )
  
  # 6. Annotate the Seurat object (if requested)
  if (!is.null(seurat_object_name)) {
    seurat_object <- get(seurat_object_name, envir = parent.frame())  # Get the Seurat object from the parent environment
    
    current_identity <- Idents(seurat_object)  # Get the current identity of the Seurat object
    current_identity_column <- NULL  # Initialize the identity column name
    
    for (col_name in colnames(seurat_object@meta.data)) {
      if (all(as.character(current_identity) == as.character(seurat_object@meta.data[[col_name]]))) {
        current_identity_column <- col_name
        break
      }
    }
    
    if (is.null(current_identity_column)) {
      stop("Error: Unable to find the current identity column in the Seurat object. ",
           "Make sure the object has an identity set using 'Idents(seurat_object)'")
    }
    
    seurat_object[[paste0("scassist_annotation_", current_identity_column)]] <- tryCatch(
      {
        plyr::mapvalues(
          x = as.character(seurat_object@meta.data[[current_identity_column]]),  # Convert to character
          from = as.character(sca_annotation$cluster_num),  # Convert to character
          to = sca_annotation$cell_type
        )
      },
      error = function(e) {
        stop("Error: An error occurred while annotating clusters. ",
             "Please check if 'sca_annotation' is available and has the correct format.", call. = FALSE)
      }
    )
    # Return the modified Seurat object
    return(seurat_object)
    #assign(seurat_object_name, seurat_object, envir = parent.frame())  # Update the Seurat object in the parent environment
  }
  else {
    # Return the annotation results if no Seurat object is provided
    return(sca_annotation)
  }
  # 7. Return the annotation results
  #return(sca_annotation)
}