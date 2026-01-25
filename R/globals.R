# Global variables for NSE (Non-Standard Evaluation)
# These are used in dplyr/data.table pipelines and ggplot2 aes()
utils::globalVariables(c(
    # tfidf.R
    "word", "var", "cluster", "n", "tf_idf",
    # plot_tfidf.R
    "level", "x.mean", "y.mean", "size", "term",
    # search_neighbors.R
    "Var1", "Var2", "similarity",
    # gpt.R - uses cluster already declared above
    NULL
))
