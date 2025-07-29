hDP_mat_joint <- function(
  n,
  group_indx,
  seen = matrix(numeric(0),ncol=4,dimnames = list(NULL,c("X", "Ntabs", "Ncusts1", "Ncusts2"))),
  c0 = 1,
  c = 1,
  P00 = runif){
    if(n==0){
      return(seen)
    }
    n_dishes <- nrow(seen)
    n_tabs <- sum(seen[,2])
    n_custs <- sum(seen[,group_indx])
    
    seen <- rbind(seen,matrix(numeric(0),nrow = n,ncol = 4))
    
    for (k in seq_len(n)) {
      if (runif(1) < c/(c + n_custs)) {
        # new table
        if (runif(1) < c0/(c0 + n_tabs)){
          # new dish
          n_dishes <- n_dishes + 1
          
          vec <- c(P00(1),1,0,0)
          vec[group_indx] <- 1
          
          seen[n_dishes,] <- vec
        }
        else{
          # old dish
          
          indx <- sample(1:n_dishes,size = 1,prob = seen[1:n_dishes,2])
          seen[indx,2] <- seen[indx,2] + 1
          seen[indx,group_indx] <- seen[indx,group_indx] + 1
        }
        n_tabs <- n_tabs + 1
      }
      else {
        # old table, old dish
        indx <- sample(1:n_dishes,size = 1,prob = seen[1:n_dishes,group_indx])
        seen[indx,group_indx] <- seen[indx,group_indx] + 1
      }
      n_custs <- n_custs + 1
    }
    return(seen[1:n_dishes,,drop=FALSE])
}