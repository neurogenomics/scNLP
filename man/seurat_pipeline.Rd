% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_pipeline.R
\name{seurat_pipeline}
\alias{seurat_pipeline}
\title{Run standardized \pkg{Seurat} pipeline}
\usage{
seurat_pipeline(...)
}
\arguments{
\item{...}{
  Arguments passed on to \code{\link[scKirby:process_seurat]{scKirby::process_seurat}}
  \describe{
    \item{\code{cluster_reduction}}{Recompute neighbors graph based on UMAP
to get clusters that best reflect UMAP space.
For this same reason, only cluster in two dimensions,
because this is the view we most often use.
That said, this may reduce the generalisability of these clusters/graph.}
  }}
}
\value{
A preprocessed \link[Seurat]{Seurat} object.
}
\description{
Run \pkg{Seurat} pipeline on \pkg{Seurat} object 
or raw \code{counts} and \code{meta.data}.\cr
Automatically performs 
\describe{
\item{\code{FindVariableFeatures}}{Variable feature selection}
\item{\code{NormalizeData}}{Data normalization}
\item{\code{RunPCA}}{PCA}
\item{\code{RunUMAP}}{UMAP}
\item{\code{FindNeighbors}}{K-nearest neighbors}
\item{\code{FindClusters}}{Clustering}
}
}
\examples{
X <- scNLP::pseudo_sce@assays@data$raw@seed[,seq_len(100)]
obj <- seurat_pipeline(obj = X)
}
