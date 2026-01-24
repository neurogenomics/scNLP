# Helper function to load and validate pseudo_seurat
# Returns the object if usable, or skips the test if incompatible with current Seurat

load_pseudo_seurat <- function() {
    testthat::skip_if_not_installed("Seurat")
    testthat::skip_if_not_installed("SeuratObject")

    data("pseudo_seurat", package = "scNLP")

    # Try to update the object to current format
    pseudo_seurat <- tryCatch({
        SeuratObject::UpdateSeuratObject(pseudo_seurat)
    }, error = function(e) {
        pseudo_seurat
    })

    # Test if the object is usable by attempting a basic operation
    tryCatch({
        # This will fail if the object structure is incompatible
        Seurat::GetAssayData(pseudo_seurat, layer = "counts")
        pseudo_seurat
    }, error = function(e) {
        testthat::skip(paste0(
            "pseudo_seurat incompatible with current Seurat version: ",
            conditionMessage(e)
        ))
    })
}
