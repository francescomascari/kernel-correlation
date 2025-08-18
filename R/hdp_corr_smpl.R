hdp_corr_smpl <- function(
    seen = matrix(numeric(0),
                  ncol = 3,
                  dimnames = list(NULL, c("X", "Ncusts1", "Ncusts2"))),
    c0 = 1,
    c = 1,
    bsln = runif,
    M = 10000,
    ...) {
  # ---------------------------------------------------------------------------
  # Compute the kernel correlation for the hierarchical DP (hDP) model
  # using the "sampling" method specified in the supplementary material
  # of the manuscript.
  # If `seen` is supplied, the computation is a posteriori,
  # conditionally on the counts in that matrix.
  # ---------------------------------------------------------------------------
  #
  # Arguments
  #   seen        : numeric matrix with columns
  #                 (using the restaurant franchise metaphor).
  #                   X        – unique dish values.
  #                   Ncusts1  – frequency of customers in group 1.
  #                   Ncusts2  – frequency of customers in group 2.
  #                 If empty (default), the computation is a priori.
  #
  #   c0, c       : numeric – concentration parameters of the hDP (default 1).
  #
  #   bsln        : function – baseline measure generator (default `runif`).
  #
  #   M           : integer – number of 2x2 samples from the model.
  #
  #   ...         : extra arguments: kernel and its parameters.
  #
  # Returns:
  #   numeric – the value of the kernel correlation.

  # store the unique values of the dishes
  X_vec <- seen[, "X"]

  # store the frequencies of the unique values for the first group
  n1_vec <- seen[, "Ncusts1"]
  n1 <- sum(n1_vec) # total number of customers for the first group

  # store the frequencies of the unique values for the second group
  n2_vec <- seen[, "Ncusts2"]
  n2 <- sum(n2_vec) # total number of customers for the second group

  # initialize by assigning one table per each customer for each group
  seen_q1_mat <- lapply(n1_vec, FUN = function(times) {rep(1, times)})
  seen_q2_mat <- lapply(n2_vec, FUN = function(times) {rep(1, times)})

  # if at least one group has some values
  if (n1 != 0 || n2 != 0) {
    # store the frequencies of unique values for the tables across both groups
    l_vec <- sapply(seen_q1_mat, length) + sapply(seen_q2_mat, length)
  } else {
    # assign an empty vector
    l_vec <- c()
  }

  # pre-allocate the vectors to store the sampled values
  X11_vals <- c()
  X12_vals <- c()
  X21_vals <- c()
  X22_vals <- c()

  for (. in seq_len(M)) {
    # we impose the table disposition obained from the Gibbs update
    seen_full <- cbind(seen, "Ntabs" = l_vec)

    # generate the 2x2 observations from the model with given `seen_full`
    joint_vals <- hdp_sampler(2, 2, seen = seen_full, smpl_method = "alt", start = 1, c0 = c0, c = c, P00 = bsln)$new_X

    # store the sampled values accordingly
    X11_vals <- c(X11_vals,joint_vals[joint_vals[, "group"] == 1, 1][1])
    X12_vals <- c(X12_vals,joint_vals[joint_vals[, "group"] == 1, 1][2])
    X21_vals <- c(X21_vals,joint_vals[joint_vals[, "group"] == 2, 1][1])
    X22_vals <- c(X22_vals,joint_vals[joint_vals[, "group"] == 2, 1][2])

    # Gibbs update
    # if at least one group has some values
    if (n1 != 0 || n2 != 0) {
      # store the frequencies of unique values for the tables across both groups
      seen_q1_mat <- gibbs_tabs(l_vec, seen_q1_mat, c0 = c0, c = c)
      l_vec <- sapply(seen_q1_mat, length) + sapply(seen_q2_mat, length)
      seen_q2_mat <- gibbs_tabs(l_vec, seen_q2_mat, c0 = c0, c = c)
      l_vec <- sapply(seen_q1_mat, length) + sapply(seen_q2_mat, length)
    } else {
      # assign an empty vector
      l_vec <- c()
    }
  }

  # compute kernel matrices
  outer_mat12 <- do_outer_mat(X11_vals, X21_vals, ...)
  outer_mat11 <- do_outer_mat(X11_vals, X12_vals, ...)
  outer_mat22 <- do_outer_mat(X21_vals, X22_vals, ...)

  # compute covariance and variances
  # NOTE: no division by (M-1) as it simplifies in the ratio
  cov <- sum(diag(outer_mat12)) - sum(outer_mat12) / M
  var1 <- sum(diag(outer_mat11)) - sum(outer_mat11) / M
  var2 <- sum(diag(outer_mat22)) - sum(outer_mat22) / M

  return(cov / sqrt(var1 * var2))
}
