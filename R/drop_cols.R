drop_cols <- function(df,
                      cols){
  for(x in cols){
    if(x %in% names(df)) df[[x]] <- NULL
  } 
  return(df)
}