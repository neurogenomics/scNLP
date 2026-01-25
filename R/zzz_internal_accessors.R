#' Get observations metadata from Seurat object
#'
#' Internal function to extract cell metadata from a Seurat object.
#' Replaces scKirby::get_obs for Bioconductor compatibility.
#'
#' @param obj A Seurat object.
#' @param verbose Print messages.
#' @returns A data.frame of cell metadata.
#' @keywords internal
#' @noRd
get_obs_internal <- function(obj, verbose = TRUE) {
    if (!methods::is(obj, "Seurat")) {
        stop("obj must be a Seurat object")
    }
    return(obj@meta.data)
}

#' Get dimensional reductions from Seurat object
#'
#' Internal function to extract dimensional reduction embeddings from a Seurat object.
#' Replaces scKirby::get_obsm for Bioconductor compatibility.
#'
#' @param obj A Seurat object.
#' @param keys Character vector of reduction names to extract (e.g., "umap", "pca").
#' @param verbose Print messages.
#' @returns A named list of embedding matrices.
#' @keywords internal
#' @noRd
get_obsm_internal <- function(obj, keys = NULL, verbose = TRUE) {
    if (!methods::is(obj, "Seurat")) {
        stop("obj must be a Seurat object")
    }

    all_keys <- SeuratObject::Reductions(obj)

    if (is.null(keys)) {
        keys <- all_keys
    } else {
        # Case-insensitive matching
        keys_lower <- tolower(keys)
        all_keys_lower <- tolower(all_keys)
        matched <- all_keys[match(keys_lower, all_keys_lower)]
        matched <- matched[!is.na(matched)]
        if (length(matched) == 0) {
            avail <- paste(all_keys, collapse = ", ")
            stop("No matching reductions found. Available: ", avail)
        }
        keys <- matched
    }

    if (verbose) {
        messager("Extracting obsm from Seurat:", paste(keys, collapse = ", "))
    }

    obsm <- lapply(stats::setNames(keys, keys), function(k) {
        as.data.frame(SeuratObject::Embeddings(obj, reduction = k))
    })

    return(obsm)
}

#' Set observations metadata in Seurat object
#'
#' Internal function to set cell metadata in a Seurat object.
#' Replaces scKirby::set_obs for Bioconductor compatibility.
#'
#' @param obj A Seurat object.
#' @param obs A data.frame of cell metadata with rownames matching cell names.
#' @param verbose Print messages.
#' @returns The modified Seurat object.
#' @keywords internal
#' @noRd
set_obs_internal <- function(obj, obs, verbose = TRUE) {
    if (!methods::is(obj, "Seurat")) {
        stop("obj must be a Seurat object")
    }

    if (verbose) {
        messager("Setting cell metadata (obs) in obj.")
    }

    # Ensure obs is a data.frame with proper rownames
    obs <- as.data.frame(obs)

    # Match to object cell names
    cell_names <- colnames(obj)
    if (!all(cell_names %in% rownames(obs))) {
        warning("Some cells in obj not found in obs. Subsetting to matching cells.")
        obs <- obs[rownames(obs) %in% cell_names, , drop = FALSE]
    }

    obj@meta.data <- obs[cell_names, , drop = FALSE]
    return(obj)
}
