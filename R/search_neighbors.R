#' Search nearest neighbors 
#' 
#' Search [shared] K-nearest neighbor graph to find the
#'  samples that are most similar to those matching a substring search.
#'  
#' @param seurat \pkg{Seurat} object. 
#' @param graph_name Name of the graph to use. 
#' If none provided, will use the last graph available.
#' If no graphs are available, new ones will be computed using \code{Seurat::FindNeighbors}. 
#' @param var1_search Substring search term to filter var1 by. 
#' If a vector is supplied instead, this will be interpreted as an "or" query. 
#' @param label_col \code{meta.data} column used to name the rows/columns of the graph.
#' \code{label_col} will also be used in the search for \code{var1_search} substring.
#' @param var2_group Substring search term to filter var2 by, 
#' according to  \code{group_col}.
#' @param group_col \code{meta.data} column used to filter var2 
#' when \code{var2_group} is used. 
#' @param max_neighbors The max number of neighbors (var2) per term (var1).   
#' @param add_original_names Add original names into the results. 
#' This can be useful when var1 names are forced to be unique internally. 
#' @param verbose Whether to print messages. 
#'      
#' @examples 
#' library(scNLP)
#' data("pseudo_seurat")
#' 
#' ### No group filter
#' top_neighbors <- search_neighbors(seurat = pseudo_seurat, 
#'                                   var1_search = "purkinje", 
#'                                   max_neighbors=5)
#' ### With group filter
#' top_neighbors2 <- search_neighbors(seurat = pseudo_seurat,
#'                                   var1_search = "purkinje",
#'                                   var2_group = "human",
#'                                   group_col = "species",
#'                                   max_neighbors=5)
#' @export                    
search_neighbors <- function(seurat,
                             graph_name=NULL,
                             var1_search=NULL,  
                             label_col=NULL,
                             var2_group=NULL,
                             group_col=NULL, 
                             max_neighbors=Inf, 
                             add_original_names=T,
                             verbose=T){ 
  if(is.null(names(seurat@graphs))){ 
    if(!"pca" %in% names(seurat@reductions)){
      if(length(Seurat::VariableFeatures(seurat))==0){
        messager("No variable features detected. Computing",v=verbose)
        seurat <- Seurat::FindVariableFeatures(seurat)
      }
      messager("No PCA detected. Computing",v=verbose)
      seurat <- Seurat::NormalizeData(seurat)
      seurat <- Seurat::ScaleData(seurat)
      seurat <- Seurat::RunPCA(seurat)
    }
    messager("No graphs detected. Computing.",v=verbose)
    seurat <- Seurat::FindNeighbors(seurat)
  }
  
  if(length(var1_search)>1) var1_search <- paste(var1_search,collapse = "|")
  if(is.null(graph_name)) graph_name <- rev(names(seurat@graphs))[1] 
  if(any(graph_name %in% names(seurat@graphs))){
    messager("Using graph:",graph_name,v=verbose)
  }
  fgraph <- seurat@graphs[[graph_name]]
  
  if(is.null(label_col)){
    sample_names <- rownames(fgraph)
  }else{
    sample_names <- make.unique(seurat@meta.data[[label_col]])
  }  
  seurat@meta.data$sample_names <- sample_names   
  names_dict <- setNames(rownames(fgraph), sample_names) 
  
  fgraph <- fgraph |> 
    `row.names<-`(sample_names) |>
    `colnames<-`(sample_names)
  
  if(!is.null(var1_search)){
    messager("+ Filtering results by `var1_search`:",var1_search,v=verbose)
    targets1 <- grep(var1_search,sample_names, ignore.case = T, value = T)
    if(length(targets1)>0){
      messager("+",length(targets1),"entries matching `var1_search` identified.",v=verbose)
      fgraph <- fgraph[targets1,]
    } else {stop("0 entries in `label_col` match the substring search for `var1_search`")} 
  }
  
  if(!is.null(group_col)){
    messager("+ Filtering results by `var2_group`:",var2_group,v=verbose)  
    targets2 <- seurat@meta.data[grepl(var2_group, seurat@meta.data[[group_col]], ignore.case = T),
                                 "sample_names"]
    if(length(targets2)>0){
      messager("+",length(targets2),"entries matching `var2_group` identified.",v=verbose)
      fgraph <- fgraph[,targets2]
    } else {stop("0 entries in `group_col` match the substring search for `var2_group`")} 
  }
  
  top_candidates <-  fgraph |>
    Matrix::as.matrix() |>
    reshape2::melt(value.name = "similarity") |>
    data.table::data.table() |>
    dplyr::mutate_at(c("Var1","Var2"), as.character) |>
    subset(Var1!=Var2) |>  
    subset(similarity>0) |>
    dplyr::group_by(Var1) |> 
    dplyr::slice_max(order_by = similarity, n = max_neighbors) |> 
    dplyr::arrange(desc(similarity)) |> 
    data.table::data.table()    
  
  if(!is.null(group_col)){
    messager("+ Adding `group_col` to results")
    keys <- setNames(seurat@meta.data[[group_col]], seurat@meta.data$sample_names) 
    top_candidates$Var2_group <- keys[top_candidates$Var2]
  } 
  if(add_original_names){
    messager("+ Adding original names to results")
    top_candidates$Var1_id <- names_dict[top_candidates$Var1]
    top_candidates$Var2_id <- names_dict[top_candidates$Var2]
  }
 
  messager("+ Returning",nrow(top_candidates),"pair-wise similarities.",v=verbose)
  return(top_candidates)
}
