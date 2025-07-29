cond_law_tabs2_log <- function(
  n1_vec,
  n2_vec,
  c0 = 10,
  c = 10){
    stopifnot(length(n1_vec)==length(n2_vec))
    k <- length(n1_vec)
    
    l1_df <- expand.grid(lapply(n1_vec,my_seq_len))
    colnames(l1_df) <- c(paste("l1_",seq_len(k)))
    l1_size <- nrow(l1_df)
    
    l2_df <- expand.grid(lapply(n2_vec,my_seq_len))
    colnames(l2_df) <- c(paste("l2_",seq_len(k)))
    l2_size <- nrow(l2_df)
    
    l_df <- data.frame(matrix(numeric(0),ncol=k))
    colnames(l_df) <- c(paste("l_",seq_len(k)))
    
    un_prob <- rep(numeric(0),l1_size*l2_size)
    for(i in seq_len(l1_size)){
      for(j in seq_len(l2_size)){
        l1_vec <- l1_df[i,]
        l2_vec <- l2_df[j,]
        
        l_df[(i-1)*l1_size + j,] <- l1_vec + l2_vec
        
        un_prob[(i-1)*l1_size + j] <- log_un_prob_tabs2(n1_vec,n2_vec,l1_vec,l2_vec,c0,c)
      }
    }
    
    l_df <- aggregate(un_prob,by=l_df,FUN=logSumExp)
    colnames(l_df)[k+1] <- "Prob"
    
    l_df$Prob <- exp(l_df$Prob - logSumExp(l_df$Prob))
    
    return(l_df)
}