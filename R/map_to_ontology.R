#### UNDER CONSTRUCTION ####
# map_to_ontology <- function(){
#   requireNamespace("celldex")
#   requireNamespace("ontoProc")
#   #### Very simple function where names must match exactly in reference, but some useful datasets ####
#   #ontoProc::  https://master.bioconductor.org/packages/release/bioc/vignettes/ontoProc/inst/doc/ontoProc.html
#   library(ontoProc)
#   hpca_map = read.csv(system.file("extdata/hpca.csv", package="ontoProc"), strings=FALSE)
#   names(hpca_map) = c("informal", "formal") 
#   library(celldex)
#   hpca_sce = celldex::HumanPrimaryCellAtlasData()
#   hpca_sce$label.ont_orig <- hpca_sce$label.ont
#   hpca_sce = ontoProc::bind_formal_tags(se = hpca_sce, 
#                                         informal = "label.fine",
#                                         tagmap = hpca_map)
#   length(unique(hpca_sce$label.ont_orig))
#   length(unique(hpca_sce$label.ont)) 
# }
# 
# 
