hDP_XT2 <- function(
  n1,
  n2,
  smpl_method = "block",
  start = 1,
  seen_XT2 = matrix(numeric(0),ncol=3,dimnames = list(NULL,c("X", "T", "group"))),
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
    
    if (smpl_method == "block") {
      seen_XT2 <- hDP_XT2_help(n_1st,group_indx = indx_1st, seen = seen_XT2, c0 = c0, c = c,P00 = P00)
      seen_XT2 <- hDP_XT2_help(n_2nd,group_indx = indx_2nd, seen = seen_XT2, c0 = c0, c = c,P00 = P00)
    }
    else if (smpl_method == "alt") {
      n_min <- min(n1,n2)
      n_max <- max(n1,n2)
      max_indx <- which.max(c(n1,n2))
      for (k in seq_len(n_min)) {
        seen_XT2 <- hDP_XT2_help(1,group_indx = indx_1st, seen = seen_XT2, c0 = c0, c = c,P00 = P00)
        seen_XT2 <- hDP_XT2_help(1,group_indx = indx_2nd, seen = seen_XT2, c0 = c0, c = c,P00 = P00)
      }
      seen_XT2 <- hDP_XT2_help(n_max - n_min,group_indx = max_indx, seen = seen_XT2, c0 = c0, c = c,P00 = P00)
    }
    return(seen_XT2)
}