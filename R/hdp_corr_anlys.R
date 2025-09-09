hdp_corr_anlys <- function(
    seen = matrix(numeric(0),
                  ncol = 3,
                  dimnames = list(NULL, c("X", "Ncusts1", "Ncusts2"))),
    c0 = 1,
    c = 1,
    R = 1000,
    bsln = runif,
    M = 10000,
    ...) {
  # ---------------------------------------------------------------------------
  # Compute the kernel correlation for the hierarchical DP (hDP) model
  # using the "analytics" method specified in the supplementary material
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
  #   R           : integer – number of runs of the Gibbs sampler to
  #                           integrate out the tables.
  #
  #   bsln        : function – baseline measure generator (default `runif`).
  #
  #   M           : integer – size of the sample from the baseline measure.
  #
  #   ...         : extra arguments: kernel and its parameters.
  #
  # Returns:
  #   numeric – the value of the kernel correlation.

  # sample M values from `bsln`
  bsln_smpl <- bsln(M)

  # store the unique values of the dishes
  X_vec <- seen[, "X"]
  K <- length(X_vec)

  # store the frequencies of the unique values for the first group
  n1_vec <- seen[, "Ncusts1"]
  n1 <- sum(n1_vec) # total number of customers for the first group

  # store the frequencies of the unique values for the second group
  n2_vec <- seen[, "Ncusts2"]
  n2 <- sum(n2_vec) # total number of customers for the second group

  # pre-compute the kernel matrices for the computation of the integrals
  outer_mat_bsln <- do_outer_mat(vals1 = bsln_smpl, vals2 = bsln_smpl, ...)
  outer_mat_dishes <- do_outer_mat(vals1 = X_vec, vals2 = X_vec, ...)
  outer_mat_bsln_dishes <- do_outer_mat(vals1 = bsln_smpl, vals2 = X_vec, ...)

  # pre-compute the kernel vectors for the computation of the integrals
  diag_vec_bsln <- diag(outer_mat_bsln)
  diag_vec_dishes <- diag(outer_mat_dishes)

  # compute the double integral of the kernel k(x,y) w.r.t.
  # the baseline measure both in x and in y
  int_cross_bsln <- sum(outer_mat_bsln) / M^2

  # compute the integral of the kernel k(x,x) w.r.t
  # the baseline measure in x
  int_diag_bsln <- sum(diag_vec_bsln) / M

  # if at least one customer in group 1
  if (n1 != 0) {

    # compute the integral of the kernel k(x,y) w.r.t.
    # the empirical distribution of the dishes
    # according to the frequency of customers in group 1 both in x and in y
    int_cross_custs1 <- quad.form(outer_mat_dishes, n1_vec) / n1^2

    # compute the integral of the kernel k(x,y) w.r.t.
    # the baseline measure in x,
    # the empirical distribution of the dishes
    # according to the frequency of customers in group 1 in y
    int_cross_bsln_custs1 <- sum(colSums(outer_mat_bsln_dishes) * n1_vec) / (M * n1)

    # compute the integral of the kernel k(x,x) w.r.t
    # the empirical distribution of the dishes
    # according to the frequency of customers in group 1 in x
    int_diag_custs1 <- sum(diag_vec_dishes * n1_vec) / n1
  } else {

    # all the integrals involving the empirical distribution of the dishes
    # according to the frequency of the customers in group 1 are null
    int_cross_custs1 <- 0
    int_cross_bsln_custs1 <- 0
    int_diag_custs1 <- 0
  }

  # third component of the variance for group 1
  V1_3 <- n1^2 / ((c + n1 + 1) * (c + n1)^2) * (int_diag_custs1 - int_cross_custs1)

  # if at least one customer in group 2
  if (n2 != 0) {

    # compute the integral of the kernel k(x,y) w.r.t.
    # the empirical distribution of the dishes
    # according to the frequency of customers in group 2 both in x and in y
    int_cross_custs2 <- quad.form(outer_mat_dishes, n2_vec) / n2^2

    # compute the integral of the kernel k(x,y) w.r.t.
    # the baseline measure in x,
    # the empirical distribution of the dishes
    # according to the frequency of customers in group 2 in y
    int_cross_bsln_custs2 <- sum(colSums(outer_mat_bsln_dishes) * n2_vec) / (M * n2)

    # compute the integral of the kernel k(x,x) w.r.t
    # the empirical distribution of the dishes
    # according to the frequency of customers in group 2 in x
    int_diag_custs2 <- sum(diag_vec_dishes * n2_vec) / n2
  } else {

    # all the integrals involving the empirical distribution of the dishes
    # according to the frequency of customers in group 2 are null
    int_cross_custs2 <- 0
    int_cross_bsln_custs2 <- 0
    int_diag_custs2 <- 0
  }

  # third component of the variance for group 2
  V2_3 <- n2^2 / ((c + n2 + 1) * (c + n2)^2) * (int_diag_custs2 - int_cross_custs2)

  # pre-allocate the values for which we have to integrate out the tables
  mean_W0 <- 0
  mean_W_vec <- rep(0, K)
  V0_1 <- 0
  V1_1 <- 0
  V1_2 <- 0
  V2_1 <- 0
  V2_2 <- 0

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

  # Gibbs iterations
  for (. in seq_len(R)) {
    # total number of tables across both groups
    l_tot <- sum(l_vec)

    # store the vector of weights W and build its mean value
    W0 <- c0 / (c0 + l_tot)
    W_vec <- l_vec / (c0 + l_tot)
    mean_W0 <- mean_W0 + W0 / R
    mean_W_vec <- mean_W_vec + W_vec / R

    # compute the integral of the kernel k(x,y) w.r.t.
    # the empirical distribution of the dishes
    # according to the frequency of tables across groups both in x and in y
    int_cross_tabs <- quad.form(outer_mat_dishes, W_vec)

    # compute the integral of the kernel k(x,y) w.r.t.
    # the baseline measure in x,
    # the empirical distribution of the dishes
    # according to the frequency of tables across both groups in y
    int_cross_bsln_tabs <- sum(colSums(outer_mat_bsln_dishes) * W_vec) / M

    # compute the integral of the kernel k(x,x) w.r.t
    # the empirical distribution of the dishes
    # according to the frequency of tables across both groups in y
    int_diag_tabs <- sum(diag_vec_dishes * W_vec)

    # first component of the posterior variance of the baseline
    V0_1 <- V0_1 + (W0 * int_diag_bsln + int_diag_tabs) / (c0 + l_tot + 1) + (1 - 1 / (c0 + l_tot + 1)) * (W0^2 * int_cross_bsln + 2 * W0 * int_cross_bsln_tabs + int_cross_tabs)

    # if at least one customer in group 1
    if (n1 != 0) {

      # compute the integral of the kernel k(x,y) w.r.t.
      # the empirical distribution of the dishes
      # according to the frequency of the tables across both groups in x,
      # the empirical distribution of the dishes
      # according to the frequency of the customers in group 1 in y
      int_cross_tabs_custs1 <- drop(quad.3form(outer_mat_dishes, W_vec, n1_vec)) / n1
    } else {

      # all the integrals involving the empirical distribution of the dishes
      # according to the frequency of the customers in group 1 are null
      int_cross_tabs_custs1 <- 0
    }

    # first and second components of the variance for group 1
    V1_1 <- V1_1 + c^2 / ((c + n1 + 1) * (c + n1)^2) * (1 - 1 / (c0 + l_tot + 1)) * ((W0 * int_diag_bsln + int_diag_tabs) - (W0^2 * int_cross_bsln + 2 * W0 * int_cross_bsln_tabs + int_cross_tabs))
    V1_2 <- V1_2 + c * n1 / ((c + n1 + 1) * (c + n1)^2) * ((W0 * int_diag_bsln + int_diag_tabs + int_diag_custs1) - 2 * (W0 * int_cross_bsln_custs1 + int_cross_tabs_custs1))

    # if at least one customer in group 2
    if (n2 != 0) {

      # compute the integral of the kernel k(x,y) w.r.t.
      # the empirical distribution of the dishes
      # according to the frequency of the tables across both groups in x,
      # the empirical distribution of the dishes
      # according to the frequency of the customers in group 2 in y
      int_cross_tabs_custs2 <- drop(quad.3form(outer_mat_dishes, W_vec, n2_vec)) / n2
    } else {

      # all the integrals involving the empirical distribution of the dishes
      # according to the frequency of the customers in group 2 are null
      int_cross_tabs_custs2 <- 0
    }

    # first and second components of the variance for group 2
    V2_1 <- V2_1 + c^2 / ((c + n2 + 1) * (c + n2)^2) * (1 - 1 / (c0 + l_tot + 1)) * ((W0 * int_diag_bsln + int_diag_tabs) - (W0^2 * int_cross_bsln + 2 * W0 * int_cross_bsln_tabs + int_cross_tabs))
    V2_2 <- V2_2 + c * n2 / ((c + n2 + 1) * (c + n2)^2) * ((W0 * int_diag_bsln + int_diag_tabs + int_diag_custs2) - 2 * (W0 * int_cross_bsln_custs2 + int_cross_tabs_custs2))

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
  # compute the integral of the kernel k(x,y) w.r.t.
  # the empirical distribution of the dishes
  # according to the mean frequency of the tables across both groups both in x and in y
  int_cross_mean_tabs <- quad.form(outer_mat_dishes, mean_W_vec)

  # compute the integral of the kernel k(x,y) w.r.t.
  # the baseline measure in x,
  # the empirical distribution of the dishes
  # according to the mean frequency of the tables across both groups in y
  int_cross_bsln_mean_tabs <- sum(colSums(outer_mat_bsln_dishes) * mean_W_vec) / M

  # second component of the variance of the baseline
  V0_2 <- mean_W0^2 * int_cross_bsln + 2 * mean_W0 * int_cross_bsln_mean_tabs + int_cross_mean_tabs

  # variances of the two groups
  var1 <- V1_1 / R + V1_2 / R + V1_3 + c^2 / (c + n1)^2 * (V0_1 / R - V0_2)
  var2 <- V2_1 / R + V2_2 / R + V2_3 + c^2 / (c + n2)^2 * (V0_1 / R - V0_2)

  # covariance
  cov <- c^2 / ((c + n1) * (c + n2)) * (V0_1 / R - V0_2)

  return(cov / sqrt(var1 * var2))
}
