hDP_mat2 <- function(
  n1,
  n2,
  smpl_method = "block",
  start = 1,
  seen = matrix(numeric(0),ncol=4,dimnames = list(NULL,c("X", "Ntabs", "Ncusts1", "Ncusts2"))),
  c = 1,
  c0 = 1,
  P00 = runif){
    n_dishes <- nrow(seen)
    n_tabs <- sum(seen[,2])
    n_custs1 <- sum(seen[,3])
    n_custs2 <- sum(seen[,4])
    
    if (start == 1) {
      indx_1st <- 3
      n_1st <- n1
      
      indx_2nd <- 4
      n_2nd <- n2
    }
    else {
      indx_1st <- 4
      n_1st <- n2
      
      indx_2nd <- 3
      n_2nd <- n1
    }
    
    if (smpl_method == "block") {
      seen <- hDP_mat_joint(n_1st,group_indx = indx_1st, seen = seen, c0 = c0, c = c,P00 = P00)
      seen <- hDP_mat_joint(n_2nd,group_indx = indx_2nd, seen = seen, c0 = c0, c = c,P00 = P00)
    }
    else if (smpl_method == "alt") {
      n_min <- min(n1,n2)
      n_max <- max(n1,n2)
      max_indx <- which.max(c(n1,n2))
      for (k in seq_len(n_min)) {
        seen <- hDP_mat_joint(1,group_indx = indx_1st, seen = seen, c0 = c0, c = c,P00 = P00)
        seen <- hDP_mat_joint(1,group_indx = indx_2nd, seen = seen, c0 = c0, c = c,P00 = P00)
      }
      seen <- hDP_mat_joint(n_max - n_min,group_indx = max_indx + 2, seen = seen, c0 = c0, c = c,P00 = P00)
    }
    return(seen)
}