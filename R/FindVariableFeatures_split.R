
#' Find variable features within subsets  
#' 
#'   Extension of \code{Seurat::FindVariableFeatures} that 
#' finds variable features within each subset of the data.
#' 
#' @param seurat Seurat object.
#' @param split.by Metadata variable to split data into subsets.
#' @param nfeatures The number of variable genes to select.
#' @param count_overall Whether \code{nfeatures} should be interpreted
#'  as the number of variable genes overall (\code{TRUE}) or per subset (\code{FALSE}).
#' @param return_nested Return a named list of variable genes per subset (\code{TRUE}) 
#' or a flattened list (\code{FALSE}).
#' @param nfeatures_max The max number of features to return. 
#' Features are prioritised according to how many variable gene subsets they appear in.
#' 
#' @export
#' @examples
#' var_features <- FindVariableFeatures_split(seurat=scNLP::pseudo_seurat, split.by = 'batch')
FindVariableFeatures_split <- function(seurat,
                                       split.by,
                                       nfeatures=4000,
                                       count_overall=FALSE,  
                                       return_nested=FALSE,
                                       nfeatures_max=2000){ 
  message('Finding variable genes per subset.')
  splits <- unique(seurat@meta.data[[split.by]]) 
  if(count_overall){
    nfeatures_per_split <- nfeatures
  } else {
    nfeatures_per_split <- round(nfeatures/batches) 
  }
  seurat_split <- Seurat::SplitObject(seurat, split.by = split.by)
  
  ### Iterate FindVariableFeatures
  var_features <- lapply(splits, function(x){
    message("Finding variable features for: ",x)
    obj <- Seurat::FindVariableFeatures(seurat_split[[x]], 
                                        nfeatures=nfeatures_per_split) 
    Seurat::VariableFeatures(obj)
  }) %>% `names<-`(splits) 
  
  #### Prioritize variable genes shared across datasets
  feature_counts <- sort(table(unlist(var_features)), decreasing = TRUE) 
  top_features <- head(names(feature_counts), nfeatures_max) 
  
  if(return_nested){
    var_features <- lapply(var_features, function(x)x[x %in% top_features]) %>% `names<-`(splits)
    return(var_features)
  } else {
    return(top_features)
  } 
}