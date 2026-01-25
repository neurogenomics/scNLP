# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

scNLP is an R package for applying Natural Language Processing (NLP) techniques to single-cell omics data. The primary use case is harmonizing non-standardized cell-type labels across datasets using TF-IDF (Term Frequency-Inverse Document Frequency) analysis.

## Build and Test Commands

```bash
# Install dependencies
Rscript -e "remotes::install_deps(dependencies = TRUE)"

# Build and check package
R CMD build .
R CMD check scNLP_*.tar.gz --as-cran

# Run all tests
Rscript -e "testthat::test_local()"

# Run a single test file
Rscript -e "testthat::test_file('tests/testthat/test-run_tfidf.R')"

# Generate documentation
Rscript -e "devtools::document()"

# Build pkgdown site
Rscript -e "pkgdown::build_site()"
```

## Architecture

### Core TF-IDF Pipeline

The package centers on a three-layer TF-IDF workflow:

1. **`tfidf()`** (`R/tfidf.R`) - Low-level function that computes TF-IDF on metadata tables using `tidytext`. Takes cluster assignments and label columns, tokenizes labels, removes stop words, and returns per-cluster enriched terms.

2. **`run_tfidf()`** (`R/run_tfidf.R`) - Mid-level wrapper that handles Seurat object I/O. Extracts metadata via internal accessors, calls `tfidf()`, and merges results back into the object's `meta.data` as `enriched_words` and `tf_idf` columns.

3. **`plot_tfidf()`** (`R/plot_tfidf.R`) - High-level visualization that calls `run_tfidf()` internally, computes cluster centers, and generates ggplot2 scatter plots with density overlays and labeled enriched terms.

### Object Handling

- **`R/zzz_internal_accessors.R`** - Internal functions (`get_obs_internal`, `get_obsm_internal`, `set_obs_internal`) for extracting/setting Seurat object metadata and dimensional reductions. These replaced external `scKirby` dependencies for Bioconductor compatibility.

- **`R/get_input_dat.R`** - Prepares input data by extracting cluster assignments and reduction coordinates, handles column name conflicts.

### Seurat Integration

- **`seurat_pipeline()`** (`R/seurat_pipeline.R`) - Standard preprocessing pipeline (NormalizeData → FindVariableFeatures → ScaleData → RunPCA → RunUMAP → FindNeighbors → FindClusters).

- **`search_neighbors()`** - Find nearest neighbors in embedding space.

- **`FindVariableFeatures_split()`** - Split-based variable feature selection.

### Experimental Features

- **`gpt()`/`run_gpt()`** - GPT-based cluster summarization (requires `gptstudio` package).
- **`wordcloud_tfidf()`/`wordcloud_heatmap()`** - Alternative visualizations.

## Key Dependencies

- **Seurat/SeuratObject** - Single-cell data structures
- **tidytext** - TF-IDF computation and text tokenization
- **ggplot2/ggrepel** - Visualization
- **data.table/dplyr** - Data manipulation

## Data

The package includes `pseudo_seurat`, a test Seurat object accessible via `data("pseudo_seurat")`. Test helper in `tests/testthat/helper-seurat.R` loads this for test fixtures.

## CI/CD

Uses `rworkflows` GitHub Action configured in `.github/workflows/rworkflows.yml`. Runs R CMD check, testthat tests, coverage via codecov, and pkgdown site builds across Ubuntu (Bioc devel), macOS, and Windows.
