% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_enrichment.R
\name{SCassist_analyze_enrichment}
\alias{SCassist_analyze_enrichment}
\title{Analyzes KEGG and GO Pathway Enrichment Results and Provides Insights from a Language Model}
\usage{
# Using Google's Gemini Pro:
SCassist_analyze_enrichment(llm_server="google", organism = "human", markers, 
         pvalue = 0.05, log2FC = 1, 
         experimental_design = "Single-cell RNA sequencing",
         temperature = 0,
         max_output_tokens = 10048,
         model_G = "gemini-1.5-flash-latest",
         model_O = "llama3",
         api_key_file = "api_keys.txt",
         output_file = "G_SCassist_summary_network_data.txt", 
         html_file = "G_SCassist_summary_network_visualization.html",
         do_kegg = TRUE, do_go = TRUE, go_ont = "BP",
         model_params = list(seed = 42, temperature = 0, num_gpu = 0))
}
\arguments{
\item{llm_server}{The LLM server to use. Options are "google" or "ollama". Default is "google".}

\item{organism}{The organism of interest. Supported organisms are:
"human", "mouse", "rat", "zebrafish", "fly", "yeast", "worm", "arabidopsis",
"chicken". Default is "human".}

\item{markers}{A data frame containing the marker genes for a cluster,
with columns for 'p_val', and 'avg_log2FC'. **Gene names should be present
as row names of the data frame.** This data frame can be the output of the `FindMarkers`
function from Seurat.}

\item{pvalue}{The p-value threshold for selecting significant genes. Default is 0.05.}

\item{log2FC}{The absolute value of the average log2 fold change 
threshold for selecting significant genes. Default is 1.}

\item{experimental_design}{A brief description of the experimental design, e.g.,
"Single-cell RNA sequencing of immune cells from patients with cancer". Default is "Single-cell RNA sequencing".}

\item{temperature}{The temperature parameter for the language model. Higher
values lead to more creative and unpredictable outputs. Default is 0.}

\item{max_output_tokens}{The maximum number of tokens the language model can generate.
Default is 10024.}

\item{model_G}{Character string specifying the Gemini model to use for
analysis. Default is "gemini-1.5-flash-latest".}

\item{model_O}{Character string specifying the Ollama model to use for
analysis. Default is "llama3".}

\item{api_key_file}{The path to the file containing the API key.}

\item{output_file}{The name of the text file to save the network data to.
Default is "L_SCassist_summary_network_data.txt".}

\item{html_file}{The name of the HTML file to save the network visualization to.
Default is "L_SCassist_summary_network_visualization.html".}

\item{do_kegg}{Whether to perform KEGG enrichment analysis. Default is TRUE.}

\item{do_go}{Whether to perform GO enrichment analysis. Default is TRUE.}

\item{go_ont}{The ontology to use for GO enrichment. Options are "BP", "CC", "MF", or "ALL". Default is "BP".}

\item{model_params}{A list of parameters to be passed to the `rollama::query` function.
This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.}
}
\value{
The function prints a response from a language model that provides insights
  on the enriched pathways, their relationships, the biological system, and potential
  key genes or targets. It also saves the enrichment results to a text file
  named "enrichment_results.txt" in the current working directory.
}
\description{
This function analyzes KEGG and/or GO pathway enrichment results obtained from a list of genes
provided by the user, and then uses a language model 
to generate insights based on the results.
}
\details{
This function was written with assistance from Google's Gemini and Meta's Llama3.
}
\author{
Vijay Nagarajan, PhD, NEI/NIH
}
\keyword{GO,}
\keyword{KEGG,}
\keyword{RNA}
\keyword{enrichment,}
\keyword{language}
\keyword{model,}
\keyword{pathway}
\keyword{sequencing}
\keyword{single-cell}
