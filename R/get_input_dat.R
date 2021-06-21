get_input_dat <- function(object=NULL,  
                          cluster_var=NULL,
                          reduction="UMAP",
                          verbose=T){
  #### Extract necessary info ####
  
  # object = list(metadata=scNLP::pseudo_seurat@meta.data, embeddings=scNLP::pseudo_seurat@reductions$umap@cell.embeddings); reduction="UMAP";
  if(is.list(object) & all(c("metadata","embeddings") %in% names(object)) ){
    printer("+ Using data from list.",v=verbose)
    metadata <- object$metadata
    embeddings <- object$embeddings
    object_type <- "list"
  } 
  # object <- scNLP::pseudo_seurat
  if(class(object) %in% c("Seurat","SeuratObject")){ 
    printer("+ Extracting data from Seurat object.",v=verbose)
    metadata <- object@meta.data 
    reduction_names <- names(object@reductions)
    reduction <- grep(reduction, reduction_names, ignore.case = T, value = T)[1]
    if(length(reduction)>0 & !is.na(reduction)){
      printer("+ Using reduction:",reduction,v=verbose) 
    } else stop("+ No matching reductions found.")
    embeddings <- object@reductions[[reduction]]@cell.embeddings 
    object_type <- "Seurat"
  }
  # object <- scNLP::pseudo_sce
  if(class(object) %in% c("SingleCellExperiment","SummarizedExperiment")){
    printer("+ Extracting data from SingleCellExperiment object.",v=verbose)
    metadata <- SummarizedExperiment::colData(object) 
    reduction_names <- SingleCellExperiment::reducedDimNames(object)
    reduction <- grep(reduction, reduction_names, ignore.case = T, value = T)[1]
    if(length(reduction)>0 & !is.na(reduction)){
      printer("+ Using reduction:",reduction,v=verbose) 
    } else stop("+ No matching reductions found.")
    embeddings <- SingleCellExperiment::reducedDim(object, reduction)
    object_type <- "SingleCellExperiment"
  } 
   
  input_dat = metadata[,!startsWith(tolower(colnames(metadata)), tolower(reduction))]
  input_dat  <- cbind(input_dat, embeddings)
  
  #### Infer cluster var #### 
  cluster_var <- if(!is.null(cluster_var)){
    if(!cluster_var %in% colnames(input_dat)){
      messager("+ cluster_var='",cluster_var,"' not found in meta.data.",v=verbose)
      infer_cluster_var(df = input_dat,v=verbose)
    }else {cluster_var}
  } else {infer_cluster_var(df = input_dat)}
  input_dat$cluster <-  input_dat[[cluster_var]]
  
  #### Create key ####
  dim_key <- setNames(colnames(embeddings),c("x","y"))
  
  ### Drop TF-IDF vars ####
  if(any(c("enriched_words","tf_idf") %in% colnames(input_dat))){
    input_dat <- input_dat %>% dplyr::select(-c(enriched_words,tf_idf))
  }
  return(list(input_dat=data.frame(input_dat),
              cluster_var=cluster_var,
              dim_key=dim_key,
              object_type=object_type))
}
