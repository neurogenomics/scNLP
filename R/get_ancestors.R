#' Get ancestors to Cell Ontology IDs
#'
#' @return The original \code{meta} object with the new \code{ancestor_col}.
#'
#' @param meta Metadata
#' @param id_col Name of the column in \code{meta} with the ontology IDs.
#' @param ontology A controlled ontology object of class \code{ontology_index}.
#' @param levels_up How many levels up the ontology hierarchy
#'  should ancestors be retrieved from.
#' @param ancestor_col Name of the column where ancestor IDs will be stored.
#' @examples
#' if(require("celldex", quietly=TRUE))
#'     hpca_sce <- celldex::HumanPrimaryCellAtlasData()
#'     meta_ancest <- get_ancestors(
#'         meta = SummarizedExperiment::colData(hpca_sce),
#'         id_col = "label.ont"
#'     )
#' }
get_ancestors <- function(meta,
                          id_col = "cell_ontology_id",
                          ontology = NULL,
                          levels_up = 1,
                          ancestor_col = paste0("ancestor", levels_up)) {
    # Get Cell Ontology if not provided
    if (is.null(ontology)) {
        if (!requireNamespace("ontoProc", quietly = TRUE)) {
            stop("Package 'ontoProc' is required for get_ancestors(). ",
                 "Please install it or provide an 'ontology' argument.")
        }
        # Use ontoProc::getOnto() which is the current API
        ontology <- tryCatch(
            ontoProc::getOnto("cellOnto"),
            error = function(e) NULL
        )
        if (is.null(ontology)) {
            stop("Could not load Cell Ontology from ontoProc. ",
                 "Please provide an 'ontology' argument directly.")
        }
    }

    requireNamespace("ontoProc")
    requireNamespace("ontologyIndex")
    get_anc <- function(cl, id_list,
                        levels_up = 1,
                        levels_down = NULL,
                        CL_only = FALSE) {
        all_ancests <- unlist(lapply(id_list, function(id) {
            if (is.na(id)) {
                return(NA)
            } else {
                ancests <- ontologyIndex::get_ancestors(cl, id)
                if (CL_only) ancests <- ancests[startsWith(ancests, "CL:")]
                if (!is.null(levels_down)) {
                    res <- ancests[min(levels_down, length(ancests))]
                } else {
                    n_ancests <- length(ancests)
                    idx <- n_ancests - min(levels_up, n_ancests - 1) - 1
                    res <- ancests[idx]
                }
                if (length(res) > 0) {
                    return(res)
                } else {
                    return(NA)
                }
            }
        }))
        n_unique <- length(unique(stats::na.omit(all_ancests)))
        message("Identified ", n_unique, " unique ancestors.")
        n_unmapped <- sum(is.na(all_ancests))
        n_total <- length(id_list)
        message(n_unmapped, " / ", n_total, " rows could not be mapped.")
        return(all_ancests)
    }

    ids <- meta[[id_col]]
    CL_dict <- ontology$name[ids]
    CL_df <- data.frame(
        CL_id = ids,
        CL_id2 = names(CL_dict),
        name = unname(CL_dict)
    )
    message("Identifying ancestors for ", length(CL_dict), " terms.")
    CL_df[[ancestor_col]] <- get_anc(
        cl = ontology,
        id_list = CL_df$CL_id2,
        levels_up = levels_up
    )
    return(CL_df)
}
