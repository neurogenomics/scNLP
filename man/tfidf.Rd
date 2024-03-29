% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tfidf.R
\name{tfidf}
\alias{tfidf}
\title{TF-IDF}
\usage{
tfidf(
  clusts,
  label_var = "dataset",
  cluster_var = "seurat_clusters",
  terms_per_cluster = 1,
  replace_regex = "[.]|[_]|[-]",
  force_new = FALSE,
  with_ties = FALSE
)
}
\arguments{
\item{clusts}{\code{data.frame}/\code{data.table} with the 
per-cell metadata and cluster assignments.}

\item{label_var}{Which cell metadata column to input to NLP analysis.}

\item{cluster_var}{Which cell metadata column to use to identify which 
cluster each cell is assigned to.}

\item{terms_per_cluster}{The maximum number of words to return per cluster.}

\item{replace_regex}{Characters by which to split \code{label_var} into terms
(i.e. tokens) for NLP enrichment analysis.}

\item{force_new}{If NLP results are already detected the metadata, 
set \code{force_new=TRUE} to replace them with new results.}

\item{with_ties}{Should ties be kept together? The default, \code{TRUE},
may return more rows than you request. Use \code{FALSE} to ignore ties,
and return the first \code{n} rows.}
}
\description{
Run tf-idf on a metadata table.
}
