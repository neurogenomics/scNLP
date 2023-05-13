# scNLP 0.1.2

## New features

* Replace subfunctions with `scKirby` funcs where possible.
* New func: `run_gpt`
  - Alternative approach to characterising each cluster (instead of TF-IDF. 
    Asks chatGPT to summarise the cluster.
* Soft-deprecate `seurat_pipeline` and move to `scKirby::process_seurat`.
  - Provide shallow wrapper in the meantime.

# scNLP 0.1.1

## New features

* Implemented `rworkflows`. 
* Added a `NEWS.md` file to track changes to the package.
