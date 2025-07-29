hDP_X2_mat_data <- function(
  n1,
  n2,
  smpl_method = "block",
  start = 1,
  seen = matrix(numeric(0),ncol=4,dimnames = list(NULL,c("X", "Ncusts1", "Ncusts2", "Ntabs"))),
  c0 = 1,
  c = 1,
  P00 = runif){
  # ---------------------------------------------------------------------------
  # Sample n1 and n2 observations from a two-group hierarchical DP (HDP) model.
  # If `seen` is supplied, draws are conditional on the counts in that matrix.
  # The augmented model with tables is used.
  # ---------------------------------------------------------------------------
  #
  # Arguments
  #   n1, n2      : integers – numbers of observations to draw from groups 1 and 2.
  #
  #   smpl_method : character – sampling schedule
  #                 * "block" – sample all from one group, then all from the other,
  #                             starting with the group indicated by `start`.
  #                 * "alt"   – alternate single observations, starting with `start`.
  #
  #   start       : integer (1 or 2) – group to draw the first observation from.
  #
  #   seen        : numeric matrix with columns,
  #                 (using the restaurant franchise metaphor)
  #                   X        – unique sampled values
  #                   Ncusts1  – frequency in group 1
  #                   Ncusts2  – frequency in group 2
  #                   Ntabs    – number of tables across the 
  #                 If empty (default), sampling is unconditional.
  #
  #   c0, c       : numeric – concentration parameters of the HDP.
  #
  #   P00         : function – baseline measure generator (default `runif`).
  #
  # Returns
  #   list(
  #     seen   = updated `seen` matrix (same column structure as above),
  #     new_X  = matrix of newly sampled rows with columns:
  #                X     – sampled value
  #                group – group index (1 or 2)
  #   )
  
  # assign the first and the second indices and the corresponding cardinality
  # based on the value of `start`.
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
  
  # initialize the matrix of newly sampled values
  new_X <- matrix(numeric(0),ncol=2,dimnames = list(NULL,c("X", "group")))
  
  if (smpl_method == "block") { # if the sampling method is "block"
    
    # sample all values from the first group
    group1_updt <- hDP_mat_X2_help(n_1st,group_indx = indx_1st, seen = seen, c0 = c0, c = c,P00 = P00)
    seen <- group1_updt$seen
    new_X <- rbind(new_X,group1_updt$new_X)
    
    # sample all values from the second group
    group2_updt <- hDP_mat_X2_help(n_2nd,group_indx = indx_2nd, seen = seen, c0 = c0, c = c,P00 = P00)
    seen <- group2_updt$seen
    new_X <- rbind(new_X,group2_updt$new_X)
  }
  else if (smpl_method == "alt") { # if the sampling method is "alt"
    
    n_min <- min(n1,n2) # minimal number of values to sample
    n_max <- max(n1,n2) # maximal number of values to sample
    max_indx <- which.max(c(n1,n2)) # group corresponding to the maximal number of values to samples
    
    
    for (k in seq_len(n_min)) { # until we have to sample from both groups
      # sample one value from the first group
      group1_updt <- hDP_mat_X2_help(1,group_indx = indx_1st, seen = seen, c0 = c0, c = c,P00 = P00)
      seen <- group1_updt$seen
      new_X <- rbind(new_X,group1_updt$new_X)
      
      # sample one value from the second group
      group2_updt <- hDP_mat_X2_help(1,group_indx = indx_2nd, seen = seen, c0 = c0, c = c,P00 = P00)
      seen <- group2_updt$seen
      new_X <- rbind(new_X,group2_updt$new_X)
    }
    
    # sample the remaining values from the group corresponding to the maximal number
    last_updt <- hDP_mat_X2_help(n_max - n_min,group_indx = max_indx, seen = seen, c0 = c0, c = c,P00 = P00)
    seen <- last_updt$seen
    new_X <- rbind(new_X,last_updt$new_X)
  }
  return(list(seen=seen,new_X=new_X))
}