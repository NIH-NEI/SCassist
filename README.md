
<img align="left" src="inst/docs/scassist-logo.png"/>
&nbsp;

----
* Authors: Vijay Nagarajan PhD
* Affiliation: Laboratory of Immunology, NEI/NIH
* Contact: nagarajanv@nih.gov
----
**SCassist** is an R package that utilizes a combination of statistical calculations and LLM-based insights to guide users through the complex process of single-cell RNA-seq data analysis. The package aims to provide recommendations, annotations, and interpretations, leading to efficient and insightful results.
<br><br>
Publication: biorxiv prepring PDF here: https://www.biorxiv.org/content/10.1101/2025.04.22.650107v2

### **Features:**

* Recommendations for quality control and filtering
* Recommendations for the most appropriate Normalization method
* Analysis of variable features to identify processes driving the cell-to-cell variation
* Analysis of top PCs to glean insights in to processes driving the system
* Recommendations for optimal number of PCs to use for downstream analysis
* Recommendations for suitable range of resolution values for clustering
* Marker gene analysis and cell type prediction with detailed reasoning
* KEGG pathway and GO enrichment analysis and integration, providing deeper insights in to system understanding

### **Benefits:**

* **Automated Recommendations:** Receive tailored recommendations for key parameters and analysis choices based on your specific data characteristics.
* **LLM-Powered Insights:** Leverage the power of LLMs to interpret complex data, uncover hidden patterns, and generate insightful summaries and explanations.
* **Intuitive Interface:** Easy-to-use functions with clear documentation and examples, making advanced analysis accessible to researchers of all levels.
* **Confidential:** Option to use a local LLM server to keep your data and analysis confidential.
* **Cost effective:** If taking advantage of google or openai models, use pay-as-you-go low cost API call options.

### **Installation:**

```R
# Install the necessary packages
install.packages("visNetwork")
install.packages("httr")
BiocManager::install("clusterProfiler")

# Install the devtools package if you don't have it
install.packages("devtools")

# Install SCassist from GitHub
devtools::install_github("NIH-NEI/SCassist")
```
**LLM Server Setup:**
* SCassist Local Ollama Server Setup
* Install Ollama following instructions here:
https://github.com/ollama/ollama
* Start ollama desktop application or from command line (ollama serve)

```
# Install rollama package to use the local ollama llm server
install.packages("rollama")

# Download the model (in R)
pull_model("llama3")
```
* SCassist Remote Google Server Setup - obtain api-key following the instructions here:
https://ai.google.dev/gemini-api/docs/api-key

* SCassist Remote OpenAI Server Setup - obtain api-key following the instructions here:
https://platform.openai.com/api-keys

### **Example Usage:**

**Download example data:** [NK, CD4+ and CD8+ T cells from LCMV infected Ifng - CTCF binding site mutant mice - GSM6625298_scRNA_LCMV_Day4_CD4_CD8_NK_WT_filtered_feature_bc_matrix.h5](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6625299&format=file&file=GSM6625299%5FscRNA%5FLCMV%5FDay4%5FCD4%5FCD8%5FNK%5FKO%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5)


```R
# Load the SCassist and Seurat packages
library(SCassist)
library(Seurat)

# Load the downloaded example file
KO <- Read10X_h5("GSM6625298_scRNA_LCMV_Day4_CD4_CD8_NK_WT_filtered_feature_bc_matrix.h5", use.names = T)

# Create seurat object
KO <- CreateSeuratObject(counts = KO[["Gene Expression"]], names.field = 2,names.delim = "\\-")

# Set api_key_file variable
api_key_file = "api_key_from_google.txt"

# Recommend quality control filters using Gemini (online)
qc_recommendations <- SCassist_analyze_quality("KO", llm_server="google", api_key_file = api_key_file)

# Recommend quality control filters using OpenAi GPT model (online)
qc_recommendations <- SCassist_analyze_quality("KO", llm_server="openai", api_key_file = api_key_file)

# Recommend quality control filters using Llama3 (local)
qc_recommendations <- SCassist_analyze_quality("KO", llm_server="ollama")

# ...and many more functions!
```

### **Tutorial Datasets:**
* [PBMCs from Birdshot Uveitis Patient](https://github.com/PulakNath/bcr-uveitis/raw/refs/heads/main/Results/cellranger/NS7R65BBTS/cellranger_output/filtered_feature_bc_matrix.h5)
* [PBMCs from Healthy Human Control](https://github.com/PulakNath/bcr-uveitis/raw/refs/heads/main/Results/cellranger/NS3R189BTS/cellranger_output/filtered_feature_bc_matrix.h5)
* [NK, CD4+ and CD8+ T cells from LCMV infected WT mice](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6625298&format=file&file=GSM6625298%5FscRNA%5FLCMV%5FDay4%5FCD4%5FCD8%5FNK%5FWT%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5)
* [NK, CD4+ and CD8+ T cells from LCMV infected Ifng - CTCF binding site mutant mice](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6625299&format=file&file=GSM6625299%5FscRNA%5FLCMV%5FDay4%5FCD4%5FCD8%5FNK%5FKO%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5)

### **Tutorials:**

Step-by-step tutorials documenting the full workflow for the example datasets are provided below:
* [Single Cell RNA-Seq Analysis of PBMC's from a healthy human Vs an Uveitis patient, using SCassist - an AI-Based Workflow Assistant](https://nih-nei.github.io/SCassist/SCassist-Tutorial-BCRUV.html)

* [Single Cell RNA-Seq Analysis of NK, CD4+ and CD8+ T cells isolated from LCMV infected WT and Ifng - CTCF binding site mutant mice, using SCassist - an AI-Based Workflow Assistant](https://nih-nei.github.io/SCassist/SCassist-Tutorial-LCMV.html)

**Old Seurat Workflows, for comparison:**

The below workflows are the original, standard workflow versions. We used these old versions to evaluate our new SCassist based workflow.
* [Single Cell RNA-Seq Analysis of PBMC's from a healthy human Vs an Uveitis patient - Old workflow](https://pulaknath.github.io/bcr-uveitis/)

* [Single Cell RNA-Seq Analysis of NK, CD4+ and CD8+ T cells isolated from LCMV infected WT and Ifng - CTCF binding site mutant mice - Old workflow](https://nih-nei.github.io/SCassist/LCMV-standard.html)

### **Documentation:**

Detailed documentation for each function, including parameters, usage, and expected outputs, is available through the `?` help function in R. For example, run ?SCassist to known about all the included functions, run ?SCassist_analyze_quality to learn about the syntax, parameters, expected inputs, defaults and outputs about the function that analyzes the quality of your single cell data and recommends filtering options.

### **Supporting Files/Scripts:**
- [SCassist Template Prompts](https://nih-nei.github.io/SCassist/SCassist-prompts.md)
- [SCassist Semantic Similarity Evaluation Script](https://nih-nei.github.io/SCassist/bert-similarity.md)

### **License:**

The license for this package can be found in the `LICENSE` file within the package directory.
