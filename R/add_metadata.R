
add_metadata <- function(object,
                         object_type,
                         new_metadata){
  if(object_type=="Seurat") {
    object@meta.data <- new_metadata
  }
  if(object_type=="SingleCellExperiment") {
    SummarizedExperiment::colData(object) <- new_metadata
  }
  if(object_type=="list") {
    object$metadata <- new_metadata
  }
  return(object)
}
