hDP_mat_X2_help <- function(
  n,
  group_indx,
  seen = matrix(numeric(0),ncol=4,dimnames = list(NULL,c("X", "Ncusts1", "Ncusts2", "Ntabs"))),
  c0 = 1,
  c = 1,
  P00 = runif){
    new_X <-  matrix(numeric(0),ncol=2,dimnames = list(NULL,c("X", "group")))
    
    if(n==0){
      return(list(seen=seen,new_X=new_X))
    }
  
    n_dishes <- nrow(seen)
    
    n_tabs <- sum(seen[,"Ntabs"])
    
    n_custs <- sum(seen[,group_indx+1])
    
    
    for (k in seq_len(n)) {
      if (runif(1) < c/(c + n_custs)) {
        # new table
        if (runif(1) < c0/(c0 + n_tabs)){
          # new dish
          dish <- P00(1)
          n_dishes <- n_dishes + 1
          
          vec <- c(dish,0,0,1)
          vec[group_indx+1] <- 1
          
          seen <- rbind(seen,vec)
        }
        else{
          # old dish
          indx <- sample(1:n_dishes,size = 1,prob = seen[1:n_dishes,"Ntabs"])
          dish <- seen[indx,"X"]
          
          seen[indx,"Ntabs"] <- seen[indx,"Ntabs"] + 1
          seen[indx,group_indx+1] <- seen[indx,group_indx+1] + 1
        }
        n_tabs <- n_tabs + 1
      }
      else {
        # old table, old dish
        indx <- sample(1:n_dishes,size = 1,prob = seen[1:n_dishes,group_indx+1])
        dish <- seen[indx,"X"]
        
        seen[indx,group_indx+1] <- seen[indx,group_indx+1] + 1
      }
      n_custs <- n_custs + 1
      new_X <- rbind(new_X,c(dish,group_indx))
    }
  
    return(list(seen = seen,new_X = new_X))
}
