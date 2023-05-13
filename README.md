scNLP
================
<img src='https://github.com/neurogenomics/scNLP/raw/main/inst/hex/hex.gif' title='Hex sticker for scNLP' height='300'><br>
[![License: MIT + file
LICENSE](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-blue.svg)](https://cran.r-project.org/web/licenses/MIT%20+%20file%20LICENSE)
[![](https://img.shields.io/badge/devel%20version-0.1.2-black.svg)](https://github.com/neurogenomics/scNLP)
[![](https://img.shields.io/github/languages/code-size/neurogenomics/scNLP.svg)](https://github.com/neurogenomics/scNLP)
[![](https://img.shields.io/github/last-commit/neurogenomics/scNLP.svg)](https://github.com/neurogenomics/scNLP/commits/main)
<br> [![R build
status](https://github.com/neurogenomics/scNLP/workflows/rworkflows/badge.svg)](https://github.com/neurogenomics/scNLP/actions)
[![](https://codecov.io/gh/neurogenomics/scNLP/branch/main/graph/badge.svg)](https://codecov.io/gh/neurogenomics/scNLP)
<br>
<a href='https://app.codecov.io/gh/neurogenomics/scNLP/tree/main' target='_blank'><img src='https://codecov.io/gh/neurogenomics/scNLP/branch/main/graphs/icicle.svg' title='Codecov icicle graph' width='200' height='50' style='vertical-align: top;'></a>  
<h4>  
Authors: <i>Brian Schilder</i>  
</h4>
<h5>  
README updated: <i>May-13-2023</i>  
</h5>

## Tools for applying natural language processing (NLP) techniques to single-cell (sc) omics data.

# Intro

When trying to re-analyze single-cell \[RNA-seq\] data that has
previously been annotated, the same cell-types are not usually labeled
in the same way (e.g. “Purkinje cells” vs. “purkinje neurons”
vs. “pkj_neurons”). This makes harmonizing data across multiple source
quite challenging. One solution is to re-annotate all cells yourself.
Alternatively, you can re-use the existing cell-type labels with natural
language processing (NLP).

Term frequency–inverse document frequency (**tf-idf**) is an NLP
technique to identify words or phrases that are enriched in one document
relative to some other larger set of documents.

In our case, our words are within the non-standardized cell labels and
our “documents” are the clusters. The goals is to find words that are
enriched in each cluster relative to all the other clusters. This can be
thought of as an NLP equivalent of finding gene markers for each
cluster.

Another use case is to identify whether certain metadata attributes
(e.g. dataset, species, brain region) are over-represented in some
clusters relative to others This is a quantitative way to assess
whether, for example, two or more datasets have successfully been
integrated (i.e. are well-“mixed”), or whether some clusters are more
representative of a particular anatomical region.

# [Documentation website](https://neurogenomics.github.io/scNLP/)

# [tf-idf vignette](https://neurogenomics.github.io/scNLP/articles/tf-idf.html)

# Quick examples

``` r
library(scNLP) 
data("pseudo_seurat")
```

## td-idf annotation

`seurat_tfidf` will run **tf-idf** on each cluster and put the results
in the **enriched_words** and **tf_idf** cols of the `meta.data`.

``` r
pseudo_seurat <- run_tfidf(obj = pseudo_seurat,
                           reduction = "UMAP",
                           cluster_var = "cluster",
                           label_var = "celltype") 
```

    ## Extracting obsm from Seurat.

    ## + Dropping 2 conflicting obs variables: UMAP.1, UMAP.2

    ## Loading required namespace: tidytext

    ## Setting cell metadata (obs) in obj.

``` r
head(pseudo_seurat@meta.data)
```

    ##                         cluster       batch species     dataset celltype label
    ## human.DRONC_human.ASC1        5 DRONC_human   human DRONC_human     ASC1  ASC1
    ## human.DRONC_human.ASC2        5 DRONC_human   human DRONC_human     ASC2  ASC2
    ## human.DRONC_human.END         9 DRONC_mouse   mouse DRONC_mouse      END   END
    ## human.DRONC_human.exCA1       0 DRONC_human   human DRONC_human    exCA1 exCA1
    ## human.DRONC_human.exCA3       0 DRONC_human   human DRONC_human    exCA3 exCA3
    ## human.DRONC_human.exDG        0 DRONC_human   human DRONC_human     exDG  exDG
    ##                         nCount_RNA nFeature_RNA RNA_snn_res.0.8 seurat_clusters
    ## human.DRONC_human.ASC1    756.6266         1693               5               5
    ## human.DRONC_human.ASC2    766.3392         1603               5               5
    ## human.DRONC_human.END     885.2824         1645               9               9
    ## human.DRONC_human.exCA1   714.6469         1677               0               0
    ## human.DRONC_human.exCA3   634.1760         1657               0               0
    ## human.DRONC_human.exDG    659.2845         1700               0               0
    ##                             UMAP_1      UMAP_2             enriched_words
    ## human.DRONC_human.ASC1  -0.4796632  0.17629431      glia; schwann; radial
    ## human.DRONC_human.ASC2  -0.6386602 -0.05231967      glia; schwann; radial
    ## human.DRONC_human.END   -7.7066403 -1.84134831 vascular; peric; pericytes
    ## human.DRONC_human.exCA1  6.2326443  1.51104526          lpn; adpn; neuron
    ## human.DRONC_human.exCA3  6.0303471  1.47096417          lpn; adpn; neuron
    ## human.DRONC_human.exDG   5.9316036  1.49563257          lpn; adpn; neuron
    ##                                                                             tf_idf
    ## human.DRONC_human.ASC1     0.198360552120631; 0.181900967132288; 0.111766521696813
    ## human.DRONC_human.ASC2     0.198360552120631; 0.181900967132288; 0.111766521696813
    ## human.DRONC_human.END                         0.528096815017439; 0.042313284392222
    ## human.DRONC_human.exCA1 0.0527542246967963; 0.0523351433907082; 0.0428030761818744
    ## human.DRONC_human.exCA3 0.0527542246967963; 0.0523351433907082; 0.0428030761818744
    ## human.DRONC_human.exDG  0.0527542246967963; 0.0523351433907082; 0.0428030761818744

## td-idf scatter plot

You can also plot the results in reduced dimensional space (e.g. UMAP).
`plot_tfidf()` will produce a list with three items.

- `data`: The processed data used to create the plot.
- `tfidf_df`: The full per-cluster TF-IDF enrichment results.
- `plot`: The `ggplot`.

### `Seurat` input

``` r
res <- plot_tfidf(obj = pseudo_seurat, 
                  label_var = "celltype", 
                  cluster_var = "cluster", 
                  show_plot = TRUE)
```

    ## Extracting obsm from Seurat.

    ## + Dropping 2 conflicting obs variables: UMAP_1, UMAP_2

    ## Setting cell metadata (obs) in obj.

    ## Warning in ggplot2::geom_point(ggplot2::aes_string(color = color_var, size =
    ## size_var, : Ignoring unknown aesthetics: label

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

# Session Info

<details>

``` r
utils::sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] scNLP_0.1.2
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rappdirs_0.3.3                SnowballC_0.7.1              
    ##   [3] scattermore_1.0               R.methodsS3_1.8.2            
    ##   [5] SeuratObject_4.1.3            tidyr_1.3.0                  
    ##   [7] ggplot2_3.4.2                 bit64_4.0.5                  
    ##   [9] knitr_1.42                    irlba_2.3.5.1                
    ##  [11] DelayedArray_0.24.0           R.utils_2.12.2               
    ##  [13] data.table_1.14.8             KEGGREST_1.38.0              
    ##  [15] RCurl_1.98-1.12               generics_0.1.3               
    ##  [17] BiocGenerics_0.44.0           cowplot_1.1.1                
    ##  [19] RSQLite_2.3.1                 RANN_2.6.1                   
    ##  [21] future_1.32.0                 tokenizers_0.3.0             
    ##  [23] bit_4.0.5                     tzdb_0.3.0                   
    ##  [25] spatstat.data_3.0-1           xml2_1.3.4                   
    ##  [27] httpuv_1.6.9                  isoband_0.2.7                
    ##  [29] SummarizedExperiment_1.28.0   orthogene_1.5.3              
    ##  [31] xfun_0.39                     hms_1.1.3                    
    ##  [33] babelgene_22.9                evaluate_0.20                
    ##  [35] promises_1.2.0.1              fansi_1.0.4                  
    ##  [37] dbplyr_2.3.2                  igraph_1.4.2                 
    ##  [39] DBI_1.1.3                     htmlwidgets_1.6.2            
    ##  [41] spatstat.geom_3.1-0           stats4_4.2.1                 
    ##  [43] purrr_1.0.1                   ellipsis_0.3.2               
    ##  [45] dplyr_1.1.2                   ggpubr_0.6.0                 
    ##  [47] backports_1.4.1               gprofiler2_0.2.1             
    ##  [49] deldir_1.0-6                  MatrixGenerics_1.10.0        
    ##  [51] vctrs_0.6.2                   SingleCellExperiment_1.20.1  
    ##  [53] Biobase_2.58.0                here_1.0.1                   
    ##  [55] ROCR_1.0-11                   abind_1.4-5                  
    ##  [57] withr_2.5.0                   cachem_1.0.8                 
    ##  [59] grr_0.9.5                     progressr_0.13.0             
    ##  [61] sctransform_0.3.5             treeio_1.23.1                
    ##  [63] goftest_1.2-3                 cluster_2.1.4                
    ##  [65] ExperimentHub_2.6.0           ape_5.7-1                    
    ##  [67] dir.expiry_1.6.0              lazyeval_0.2.2               
    ##  [69] crayon_1.5.2                  basilisk.utils_1.10.0        
    ##  [71] crul_1.3                      spatstat.explore_3.1-0       
    ##  [73] labeling_0.4.2                pkgconfig_2.0.3              
    ##  [75] GenomeInfoDb_1.34.9           nlme_3.1-162                 
    ##  [77] pals_1.7                      badger_0.2.3                 
    ##  [79] ewceData_1.7.1                rlang_1.1.1                  
    ##  [81] globals_0.16.2                lifecycle_1.0.3              
    ##  [83] miniUI_0.1.1.1                filelock_1.0.2               
    ##  [85] httpcode_0.3.0                BiocFileCache_2.6.1          
    ##  [87] dichromat_2.0-0.1             tidytext_0.4.1               
    ##  [89] AnnotationHub_3.6.0           rprojroot_2.0.3              
    ##  [91] polyclip_1.10-4               matrixStats_0.63.0           
    ##  [93] lmtest_0.9-40                 graph_1.76.0                 
    ##  [95] Matrix_1.5-4                  aplot_0.1.10                 
    ##  [97] carData_3.0-5                 zoo_1.8-12                   
    ##  [99] whisker_0.4.1                 ggridges_0.5.4               
    ## [101] png_0.1-8                     viridisLite_0.4.2            
    ## [103] bitops_1.0-7                  R.oo_1.25.0                  
    ## [105] KernSmooth_2.23-21            Biostrings_2.66.0            
    ## [107] blob_1.2.4                    stringr_1.5.0                
    ## [109] scKirby_0.1.2                 parallelly_1.35.0            
    ## [111] spatstat.random_3.1-4         readr_2.1.4                  
    ## [113] rstatix_0.7.2                 gridGraphics_0.5-1           
    ## [115] S4Vectors_0.36.2              echodata_0.99.16             
    ## [117] ggsignif_0.6.4                scales_1.2.1                 
    ## [119] memoise_2.0.1                 magrittr_2.0.3               
    ## [121] plyr_1.8.8                    ica_1.0-3                    
    ## [123] zlibbioc_1.44.0               compiler_4.2.1               
    ## [125] echoconda_0.99.9              RColorBrewer_1.1-3           
    ## [127] fitdistrplus_1.1-11           homologene_1.4.68.19.3.27    
    ## [129] cli_3.6.1                     XVector_0.38.0               
    ## [131] listenv_0.9.0                 janeaustenr_1.0.0            
    ## [133] patchwork_1.1.2               pbapply_1.7-0                
    ## [135] MASS_7.3-60                   tidyselect_1.2.0             
    ## [137] stringi_1.7.12                highr_0.10                   
    ## [139] yaml_2.3.7                    ggrepel_0.9.3                
    ## [141] biocViews_1.66.3              grid_4.2.1                   
    ## [143] tools_4.2.1                   future.apply_1.10.0          
    ## [145] parallel_4.2.1                rworkflows_0.99.9            
    ## [147] rstudioapi_0.14               RNOmni_1.0.1                 
    ## [149] piggyback_0.1.4               gridExtra_2.3                
    ## [151] farver_2.1.1                  Rtsne_0.16                   
    ## [153] HGNChelper_0.8.1              digest_0.6.31                
    ## [155] rvcheck_0.2.1                 BiocManager_1.30.20          
    ## [157] shiny_1.7.4                   Rcpp_1.0.10                  
    ## [159] GenomicRanges_1.50.2          car_3.1-2                    
    ## [161] broom_1.0.4                   BiocVersion_3.16.0           
    ## [163] later_1.3.1                   BiocPkgTools_1.16.1          
    ## [165] RcppAnnoy_0.0.20              AnnotationDbi_1.60.2         
    ## [167] httr_1.4.5                    fauxpas_0.5.2                
    ## [169] colorspace_2.1-0              rvest_1.0.3                  
    ## [171] XML_3.99-0.14                 tensor_1.5                   
    ## [173] reticulate_1.28               rorcid_0.7.0                 
    ## [175] IRanges_2.32.0                splines_4.2.1                
    ## [177] uwot_0.1.14                   yulab.utils_0.0.6            
    ## [179] RBGL_1.74.0                   tidytree_0.4.2               
    ## [181] spatstat.utils_3.0-2          gh_1.4.0                     
    ## [183] sp_1.6-0                      basilisk_1.10.2              
    ## [185] mapproj_1.2.11                renv_0.17.3                  
    ## [187] ggplotify_0.1.0               plotly_4.10.1                
    ## [189] xtable_1.8-4                  jsonlite_1.8.4               
    ## [191] ggtree_3.6.2                  sceasy_0.0.7                 
    ## [193] ggfun_0.0.9                   R6_2.5.1                     
    ## [195] RUnit_0.4.32                  EWCE_1.9.0                   
    ## [197] pillar_1.9.0                  htmltools_0.5.5              
    ## [199] mime_0.12                     glue_1.6.2                   
    ## [201] fastmap_1.1.1                 DT_0.27                      
    ## [203] interactiveDisplayBase_1.36.0 codetools_0.2-19             
    ## [205] maps_3.4.1                    utf8_1.2.3                   
    ## [207] lattice_0.21-8                spatstat.sparse_3.0-1        
    ## [209] tibble_3.2.1                  dlstats_0.1.6                
    ## [211] curl_5.0.0                    leiden_0.4.3                 
    ## [213] zip_2.3.0                     openxlsx_4.2.5.2             
    ## [215] limma_3.54.2                  survival_3.5-5               
    ## [217] rmarkdown_2.21                desc_1.4.2                   
    ## [219] munsell_0.5.0                 GenomeInfoDbData_1.2.9       
    ## [221] reshape2_1.4.4                gtable_0.3.3                 
    ## [223] Seurat_4.3.0

</details>
