
run_tfidf <- function(object=NULL, 
                      reduction="UMAP",
                      label_var="label",
                      cluster_var="seurat_clusters",
                      replace_regex = " ",
                      terms_per_cluster=3,
                      force_new=F,
                      return_all_results=F,
                      verbose=T){
  # label_var="label";cluster_var="seurat_clusters";replace_regex = " ";terms_per_cluster=3;force_new=F
  
  #### Prepare input_dat ####
  gid <-  get_input_dat(object=object,  
                        cluster_var=cluster_var,
                        reduction=reduction,
                        verbose=verbose)
  input_dat <- gid$input_dat; dim_key <- gid$dim_key; cluster_var <- gid$cluster_var; object_type <- gid$object_type
  
  if(any(c("enriched_words","tf_idf") %in% colnames(input_dat))){
    if(force_new){
      input_dat <- input_dat %>% dplyr::select(-c(enriched_words,tf_idf))
    }else {
      message("Previous TF-IDF results detected. Use force_new=T to re-run.")
      return(seurat)
    }
  }
  tfidf_df <- tfidf(clusts = input_dat,
                    col_name =  label_var,
                    cluster_var = cluster_var,
                    replace_regex = replace_regex,
                    terms_per_cluster = terms_per_cluster)
  new_metadata <- merge(x=input_dat,
                        y=tfidf_df %>%
                          dplyr::group_by(cluster) %>%
                          dplyr::summarise(enriched_words=paste(unique(word),collapse="; "),
                                           tf_idf=paste(unique(tf_idf),collapse="; ")),
                        all.x = T,
                        by.x=cluster_var, by.y = "cluster", sort=F)
  # Make sure rows are in the right order
  new_metadata <- new_metadata[match(input_dat[[label_var]], new_metadata[[label_var]]),]
  if(sum(new_metadata[[label_var]]!=input_dat[[label_var]], na.rm = T)>0){
    stop("Seurat sample names and tfidf samples names are not aligned!")
  }
  new_metadata <- data.frame(new_metadata, row.names = row.names(input_dat))   
  new_object <- add_metadata(object = object,
                             object_type = object_type, 
                             new_metadata = new_metadata)
  if(return_all_results) {list(object=new_object, 
                               tfidf_df=tfidf_df, 
                               old_metadata=input_dat,
                               new_metadata=new_metadata, 
                               dim_key=dim_key,
                               cluster_var=cluster_var,
                               object_type=object_type)
    } else { return(new_object) }
}



