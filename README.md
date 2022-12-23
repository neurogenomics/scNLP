scNLP
================
[![](https://img.shields.io/badge/devel%20version-0.1.1-black.svg)](https://github.com/bschilder/scNLP)<br><br>
[![R build
status](https://github.com/bschilder/scNLP/workflows/rworkflows/badge.svg)](https://github.com/bschilder/scNLP/actions)
[![](https://img.shields.io/github/last-commit/bschilder/scNLP.svg)](https://github.com/bschilder/scNLP/commits/master)
[![](https://img.shields.io/github/languages/code-size/bschilder/scNLP.svg)](https://github.com/bschilder/scNLP)
[![](https://app.codecov.io/gh/bschilder/scNLP/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bschilder/scNLP)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
¶ <h4> ¶ Authors: <i>Brian Schilder</i> ¶ </h4>
<h5> ¶ README updated: <i>Dec-22-2022</i> ¶ </h5>

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
pseudo_seurat <- run_tfidf(object = pseudo_seurat,
                           reduction = "UMAP",
                           cluster_var = "cluster",
                           label_var = "celltype") 
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## [1] "+ Extracting data from Seurat object."
    ## [1] "+ Using reduction: umap"
    ## [1] "+ Dropping 2 conflicting metadata variables: UMAP.1, UMAP.2"

    ## Joining, by = "word"

    ## Joining, by = "cluster"
    ## Joining, by = "cluster"

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
res <- plot_tfidf(object = pseudo_seurat, 
                  label_var = "celltype", 
                  cluster_var = "cluster", 
                  show_plot = TRUE)
```

    ## [1] "+ Extracting data from Seurat object."
    ## [1] "+ Using reduction: umap"
    ## [1] "+ Dropping 2 conflicting metadata variables: UMAP_1, UMAP_2"

    ## Joining, by = "word"
    ## Joining, by = "cluster"
    ## Joining, by = "cluster"

    ## Warning in geom_point(aes_string(color = color_var, size = size_var, label =
    ## label_var, : Ignoring unknown aesthetics: label

    ## Warning: The dot-dot notation (`..level..`) was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `after_stat(level)` instead.
    ## ℹ The deprecated feature was likely used in the scNLP package.
    ##   Please report the issue at <]8;;https://github.com/bschilder/scNLP/issueshttps://github.com/bschilder/scNLP/issues]8;;>.

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
    ## [1] ggplot2_3.4.0      tidytext_0.4.0     SeuratObject_4.1.3 sp_1.5-1          
    ## [5] scNLP_0.1.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.4          maps_3.4.1          jsonlite_1.8.4     
    ##  [4] here_1.0.1          assertthat_0.2.1    expm_0.999-6       
    ##  [7] highr_0.10          BiocManager_1.30.19 rvcheck_0.2.1      
    ## [10] gld_2.6.6           lmom_2.9            yulab.utils_0.0.6  
    ## [13] cellranger_1.1.0    ggrepel_0.9.2       yaml_2.3.6         
    ## [16] globals_0.16.2      pillar_1.8.1        lattice_0.20-45    
    ## [19] glue_1.6.2          digest_0.6.31       RColorBrewer_1.1-3 
    ## [22] colorspace_2.0-3    htmltools_0.5.4     Matrix_1.5-3       
    ## [25] rworkflows_0.99.3   pkgconfig_2.0.3     listenv_0.9.0      
    ## [28] mvtnorm_1.1-3       scales_1.2.1        rootSolve_1.8.2.3  
    ## [31] tibble_3.1.8        proxy_0.4-27        generics_0.1.3     
    ## [34] farver_2.1.1        withr_2.5.0         cli_3.5.0          
    ## [37] magrittr_2.0.3      readxl_1.4.1        evaluate_0.19      
    ## [40] tokenizers_0.2.3    janeaustenr_1.0.0   badger_0.2.2       
    ## [43] future_1.30.0       fansi_1.0.3         parallelly_1.33.0  
    ## [46] MASS_7.3-58.1       SnowballC_0.7.0     class_7.3-20       
    ## [49] progressr_0.12.0    tools_4.2.1         data.table_1.14.6  
    ## [52] lifecycle_1.0.3     stringr_1.5.0       Exact_3.2          
    ## [55] munsell_0.5.0       isoband_0.2.7       compiler_4.2.1     
    ## [58] e1071_1.7-12        rlang_1.0.6         grid_4.2.1         
    ## [61] dichromat_2.0-0.1   rstudioapi_0.14     labeling_0.4.2     
    ## [64] rmarkdown_2.19      boot_1.3-28.1       DescTools_0.99.47  
    ## [67] gtable_0.3.1        codetools_0.2-18    DBI_1.1.3          
    ## [70] R6_2.5.1            knitr_1.41          dplyr_1.0.10       
    ## [73] fastmap_1.1.0       future.apply_1.10.0 utf8_1.2.2         
    ## [76] rprojroot_2.0.3     pals_1.7            dlstats_0.1.6      
    ## [79] desc_1.4.2          stringi_1.7.8       parallel_4.2.1     
    ## [82] Rcpp_1.0.9          mapproj_1.2.9       vctrs_0.5.1        
    ## [85] tidyselect_1.2.0    xfun_0.36

</details>
