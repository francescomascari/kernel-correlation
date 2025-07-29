mean_fun_ell <- function(
  l_df,
  FUN,
  ...){
    l_size <- nrow(l_df)
    n_col <- ncol(l_df)
    val <- 0
    for(i in seq_len(l_size)){
      l_vec <- as.numeric(l_df[i,1:n_col-1])
      val <- val + FUN(l_vec,...)*l_df$Prob[i]
    }
    return(val)
}