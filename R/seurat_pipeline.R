#' Run standardized Seurat pipeline
#'
#' Run a standard Seurat preprocessing pipeline on a Seurat object
#' or raw counts matrix.
#' Automatically performs:
#' \describe{
#' \item{\code{NormalizeData}}{Data normalization}
#' \item{\code{FindVariableFeatures}}{Variable feature selection}
#' \item{\code{ScaleData}}{Data scaling}
#' \item{\code{RunPCA}}{PCA}
#' \item{\code{RunUMAP}}{UMAP}
#' \item{\code{FindNeighbors}}{K-nearest neighbors}
#' \item{\code{FindClusters}}{Clustering}
#' }
#'
#' @param obj A Seurat object or a counts matrix.
#' @param dims Dimensions to use for UMAP and neighbors (default: seq_len(30)).
#' @param resolution Clustering resolution (default: 0.8).
#' @param verbose Print progress messages.
#' @param ... Additional arguments passed to Seurat functions.
#'
#' @returns A preprocessed \link[SeuratObject]{Seurat} object with PCA, UMAP,
#' neighbors, and cluster assignments.
#' @export
#' @import Seurat
#' @examples
#' data("pseudo_seurat")
#' # Re-run pipeline on existing object
#' obj <- seurat_pipeline(obj = pseudo_seurat)
seurat_pipeline <- function(obj,
                            dims = seq_len(30),
                            resolution = 0.8,
                            verbose = TRUE,
                            ...) {
    # Convert matrix to Seurat if needed
    if (methods::is(obj, "matrix") || methods::is(obj, "dgCMatrix")) {
        if (verbose) messager("Creating Seurat object from matrix.")
        obj <- SeuratObject::CreateSeuratObject(counts = obj, ...)
    }

    if (!methods::is(obj, "Seurat")) {
        stop("obj must be a Seurat object or a counts matrix")
    }

    # Standard pipeline
    if (verbose) messager("Running NormalizeData...")
    obj <- Seurat::NormalizeData(obj, verbose = verbose)

    if (verbose) messager("Running FindVariableFeatures...")
    obj <- Seurat::FindVariableFeatures(obj, verbose = verbose)

    if (verbose) messager("Running ScaleData...")
    obj <- Seurat::ScaleData(obj, verbose = verbose)

    if (verbose) messager("Running PCA...")
    obj <- Seurat::RunPCA(obj, verbose = verbose)

    # Adjust dims if too many requested
    max_dims <- ncol(obj@reductions$pca@cell.embeddings)
    if (max(dims) > max_dims) {
        dims <- seq_len(min(max_dims, max(dims)))
        if (verbose) messager("Adjusted dims to 1:", max(dims))
    }

    if (verbose) messager("Running UMAP...")
    obj <- Seurat::RunUMAP(obj, dims = dims, verbose = verbose)

    if (verbose) messager("Running FindNeighbors...")
    obj <- Seurat::FindNeighbors(obj, dims = dims, verbose = verbose)

    if (verbose) messager("Running FindClusters...")
    obj <- Seurat::FindClusters(obj, resolution = resolution, verbose = verbose)

    return(obj)
}
