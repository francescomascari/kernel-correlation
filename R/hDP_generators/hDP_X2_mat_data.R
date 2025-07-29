hDP_X2_mat_data <- function(
  n1,
  n2,
  smpl_method = "block",
  start = 1,
  seen = matrix(numeric(0),ncol=4,dimnames = list(NULL,c("X", "Ncusts1", "Ncusts2", "Ntabs"))),
  c = 1,
  c0 = 1,
  P00 = runif){
    if (start == 1) {
      indx_1st <- 1
      n_1st <- n1
      
      indx_2nd <- 2
      n_2nd <- n2
    }
    else {
      indx_1st <- 2
      n_1st <- n2
      
      indx_2nd <- 1
      n_2nd <- n1
    }
  
    new_X <- matrix(numeric(0),ncol=2,dimnames = list(NULL,c("X", "group")))
    
    if (smpl_method == "block") {
      group1_updt <- hDP_mat_X2_help(n_1st,group_indx = indx_1st, seen = seen, c0 = c0, c = c,P00 = P00)
      seen <- group1_updt$seen
      new_X <- rbind(new_X,group1_updt$new_X)
      
      group2_updt <- hDP_mat_X2_help(n_2nd,group_indx = indx_2nd, seen = seen, c0 = c0, c = c,P00 = P00)
      seen <- group2_updt$seen
      new_X <- rbind(new_X,group2_updt$new_X)
    }
    else if (smpl_method == "alt") {
      n_min <- min(n1,n2)
      n_max <- max(n1,n2)
      max_indx <- which.max(c(n1,n2))
      for (k in seq_len(n_min)) {
        group1_updt <- hDP_mat_X2_help(1,group_indx = indx_1st, seen = seen, c0 = c0, c = c,P00 = P00)
        seen <- group1_updt$seen
        new_X <- rbind(new_X,group1_updt$new_X)
        
        group2_updt <- hDP_mat_X2_help(1,group_indx = indx_2nd, seen = seen, c0 = c0, c = c,P00 = P00)
        seen <- group2_updt$seen
        new_X <- rbind(new_X,group2_updt$new_X)
      }
      last_updt <- hDP_mat_X2_help(n_max - n_min,group_indx = max_indx, seen = seen, c0 = c0, c = c,P00 = P00)
      seen <- last_updt$seen
      new_X <- rbind(new_X,last_updt$new_X)
    }
    return(list(seen=seen,new_X=new_X))
}
