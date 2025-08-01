hDP_XT2_help <- function(
  n,
  group_indx,
  seen_XT2 = matrix(numeric(0),ncol=3,dimnames = list(NULL,c("X", "T","group"))),
  c0 = 1,
  c = 1,
  P00 = runif){
    if(n==0){
      return(seen_XT2)
    }
    row_seen <- nrow(seen_XT2)
    
    unique_seen_XT2 <- matrix(seen_XT2[!duplicated(seen_XT2[,2]),1:2],ncol=2,dimnames=list(NULL,c("X", "T")))
    n_tabs <- length(unique_seen_XT2)
    
    n_custs <- sum(seen_XT2[,3] == group_indx)
    
    for (k in seq_len(n)) {
      if (runif(1) < c/(c + n_custs)) {
        # new table
        tab <- runif(1)
        if (runif(1) < c0/(c0 + n_tabs)){
          # new dish
          dish <- P00(1)                   
        }
        else{
          # old dish
          dish <- sample(unique_seen_XT2[,1],size = 1)
        }
        n_tabs <- n_tabs + 1
        
        seen_XT2 <- rbind(seen_XT2,c(dish,tab,group_indx))
        unique_seen_XT2 <- rbind(unique_seen_XT2,c(dish,tab))
      }
      else {
        # old table, old dish
        indx <- sample(seq_len(row_seen+k-1),size = 1,prob= as.integer(seen_XT2[,3] == group_indx))
        seen_XT2 <- rbind(seen_XT2,seen_XT2[indx,])
      }
      n_custs <- n_custs + 1
    }
    return(seen_XT2)
}