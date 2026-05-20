# SCassist development instructions for Codex

## Project type

SCassist is an R package for AI-assisted single-cell analysis workflows.

Follow the existing package structure and style:
- R package code lives in `R/`
- Tests live in `tests/testthat/`
- Documentation is roxygen-first
- User-facing functions follow the `SCassist_*` naming convention
- Existing analysis files use names such as `analyze_quality.R`, `analyze_pcs.R`, and `analyze_enrichment.R`

Do not create a Python framework or an `agents/` directory unless explicitly requested.

## Local R environment

This project uses a conda environment named `scassist`.

Do not assume `Rscript` is directly on PATH.

Use this command pattern for R commands:

/Users/puthranar/miniforge3/bin/conda run -n scassist Rscript -e '...'

## Recommended checks

Run InteractionAgent tests with:

/Users/puthranar/miniforge3/bin/conda run -n scassist Rscript -e 'devtools::test(filter = "analyze_interactions")'

Run the full test suite with:

/Users/puthranar/miniforge3/bin/conda run -n scassist Rscript -e 'devtools::/Users/puthran pa/Users/puthranar/miniforge3/bin/conda run -n scassist Rscript -e 'devscass/Users/puthranar/miniforge3/bin/conda run -n sPac/Users/puthranar/miniforge3/'


Users/puthranar/miniforge3/bin/ngUsers/puthranar/miniforge3/b refaUsers/puthranar/miniforge3/bin/ngUsers/puthranar/miniforge3/b reon style and workflow.
- Keep CellChat optional; SCassist should still load if CellChat is missing.
- Prefer bas- Prefer bas- Prefer bas- Prefer bas- Prefer bas- Prefer bas- Prefer bas- Prefer baonal d- Prefer bas- Prefer cefully if packages are unavailable.
