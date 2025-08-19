hdp_XT2mat <- function(seen_XT) {
  # ---------------------------------------------------------------------------
  # Convert the `seen_XT` matrix into the `seen_mat` matrix,
  # storing all the unique dish with the frequencies of tables and customers.
  # ---------------------------------------------------------------------------
  #
  # Arguments
  #   seen_XT     : numeric matrix with columns,
  #                 (using the restaurant franchise metaphor)
  #                   X      – label of the dish
  #                   T      – label of the table
  #                   group  – group of the observation (1/2)
  #
  # Returns
  #   numeric matrix with columns,
  #                 (using the restaurant franchise metaphor)
  #                   X        – unique sampled values
  #                   Ncusts1  – frequency in group 1
  #                   Ncusts2  – frequency in group 2
  #                   Ntabs    – number of tables across restaurants

  # store the unique dishes
  X <- unique(seen_XT[, "X"])

  # store the number of tables per dish
  Ntabs <- sapply(X, function(val) {length(unique(seen_XT[seen_XT[, "X"] == val, "T"]))})

  # store the number of customers per dish per restaurant
  Ncusts1 <- sapply(X, function(val) {sum(seen_XT[seen_XT[,"group"] == 1, "X"] == val)})
  Ncusts2 <- sapply(X, function(val) {sum(seen_XT[seen_XT[,"group"] == 2, "X"] == val)})

  return(data.frame(X = X, Ntabs = Ntabs, Ncusts1 = Ncusts1, Ncusts2 = Ncusts2))
}