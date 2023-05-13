#' Get ancestors to Cell Ontology IDs
#' 
#' @return The original \code{meta} object with the new \code{ancestor_col}.  
#' 
#' @param meta Metadata
#' @param id_col Name of the column in \code{meta} with the ontology IDs.
#' @param ontology A controlled ontology object of class \code{ontology_index}. 
#' @param levels_up How many levels up the ontology hierarchy
#'  should ancestors be retrieved from.
#' @param ancestor_col The name of the column where the ancestor IDs will be stored. 
#' @examples
#' \dontrun{
#'   hpca_sce = celldex::HumanPrimaryCellAtlasData()
#'   meta_ancest <- get_ancestors(meta=hpca_sce@colData , id_col="label.ont")
#' }
get_ancestors <- function(meta,
                          id_col="cell_ontology_id",
                          ontology = ontoProc::getCellOnto(),
                          levels_up=1,
                          ancestor_col=paste0("ancestor",levels_up)
                          ){   
  
  requireNamespace("ontoProc")
  requireNamespace("ontologyIndex")
  get_anc <- function(cl, id_list, 
                      levels_up=1, 
                      levels_down=NULL,
                      CL_only=FALSE){
    all_ancests <- unlist(lapply(id_list, function(id){
      if(is.na(id)){
        return(NA)
      }else {
        ancests <- ontologyIndex::get_ancestors(cl, id)
        if(CL_only) ancests <- ancests[startsWith(ancests, "CL:")] 
        if(!is.null(levels_down)){
          res <- ancests[min(levels_down, length(ancests))] 
        } else {
          res <- ancests[length(ancests)-min(levels_up, length(ancests)-1)-1] 
        } 
        if(length(res)>0){
          return(res)
        }else {return(NA)}
      }
     
    }))
    message("Identified ",length(unique(na.omit(all_ancests)))," unique ancestors.")
    message(sum(is.na(all_ancests))," / ",length(id_list)," rows could not be mapped.")
    return(all_ancests)
  }
  
  ids <- meta[[id_col]]
  CL_dict <- ontology$name[ids]  
  CL_df = data.frame(CL_id=ids,
                     CL_id2=names(CL_dict),
                     name=unname(CL_dict))  
  message("Identifying ancestors for ",length(CL_dict)," terms.")
  CL_df[[ancestor_col]] <- get_anc(cl = ontology, id_list = CL_df$CL_id2, levels_up = levels_up)
  return(CL_df)
}



# get_ancest_df <- function(){
#   
#   library(ontologyIndex)
#   data("hpo")
#   
#   if(!exists("ancest_df")){ 
#     ancest_1 <- lapply(unique(HPO$HPOid), function(x){
#       # print(x)
#       tryCatch(expr = {
#         terms <- get_term_property(ontology=hpo, property="ancestors", term=x, as_names=T)
#         return(terms[2]) 
#       }, 
#       error = function(e)return(NULL))
#     }) |> `names<-`(unique(HPO$HPOid))
#     ancest_2 <- lapply(unique(HPO$HPOid), function(x){
#       # print(x)
#       tryCatch(expr = {
#         terms <- get_term_property(ontology=hpo, property="ancestors", term=x, as_names=T)
#         return(terms[3]) 
#       }, 
#       error = function(e)return(NULL))
#     }) |> `names<-`(unique(HPO$HPOid))
#     
#     
#     ancest_df <- rbind(data.frame(ancestor_label=unlist(ancest_1), ancestor_lvl=1), 
#                        data.frame(ancestor_label=unlist(ancest_2),  ancestor_lvl=2)) |> 
#       tibble::rownames_to_column(var = "id") |>
#       tidyr::separate(col = "id", sep = "[.]", into=c("HPOid", "ancestor_id"))  |>
#       data.table::data.table() |>
#       data.table::melt.data.table(id.vars = c("HPOid","ancestor_lvl"), 
#                                   measure.vars = c("ancestor_id","ancestor_label")) |>
#       dplyr::mutate(variable=paste0(variable,ancestor_lvl)) |>
#       data.table::dcast.data.table(formula = HPOid ~ variable, value.var = "value")
#     ancest_df
#     
#     n_distinct(ancest_df$HPOid)  
#   }
# }





