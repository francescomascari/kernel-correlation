hdp_XT_sampler <- function(
    n1,
    n2,
    smpl_method = "block",
    start = 1,
    seen_XT = matrix(numeric(0),
                     ncol = 3,
                     dimnames = list(NULL, c("X", "T", "group"))),
    c0 = 1,
    c = 1,
    P00 = runif) {
  # ---------------------------------------------------------------------------
  # Sample n1 and n2 observations from a two-group hierarchical DP (HDP) model
  # with both the dish and the table labels.
  # If `seen_XT` is supplied, draws are conditional on the counts in that matrix.
  # The augmented model with tables is used.
  # ---------------------------------------------------------------------------
  #
  # Arguments
  #   n1, n2      : integers – numbers of observations to draw
  #                            from groups 1 and 2, respectively
  #
  #   smpl_method : character – sampling schedule
  #                 * "block" – sample all values from one group,
  #                             then all values from the other group,
  #                             starting with the group indicated by `start`.
  #                 * "alt"   – sample single observations,
  #                             alternating between groups,
  #                             starting with `start`.
  #
  #   start       : integer (1 or 2) – group to draw the first observation from.
  #
  #   seen_XT     : numeric matrix with columns,
  #                 (using the restaurant franchise metaphor)
  #                   X      – label of the dish
  #                   T      – label of the table
  #                   group  – group of the observation (1/2)
  #                 If empty (default), sampling is unconditional.
  #
  #   c0, c       : numeric – concentration parameters of the HDP.
  #
  #   P00         : function – baseline measure generator (default `runif`).
  #
  # Returns
  #   numeric matrix with columns : updated `seen_XT` matrix
  #                                 (same column structure as above).

  # assign the first and the second indices and the corresponding cardinality
  # based on the value of `start`.
  if (start == 1) {
    indx_1st <- 1
    n_1st <- n1

    indx_2nd <- 2
    n_2nd <- n2
  } else {
    indx_1st <- 2
    n_1st <- n2

    indx_2nd <- 1
    n_2nd <- n1
  }

  if (smpl_method == "block") { # if the sampling method is "block"
    # sample all values from the first group
    seen_XT <- hdp_XT_sampler_help(n_1st, group_indx = indx_1st, seen_XT = seen_XT, c0 = c0, c = c, P00 = P00)

    # sample all values from the second group
    seen_XT <- hdp_XT_sampler_help(n_2nd, group_indx = indx_2nd, seen_XT = seen_XT, c0 = c0, c = c, P00 = P00)
  } else if (smpl_method == "alt") { # if the sampling method is "alt"

    n_min <- min(n1, n2) # minimal number of values to sample
    n_max <- max(n1, n2) # maximal number of values to sample
    max_indx <- which.max(c(n1, n2)) # group corresponding to the maximal number of values to sample
    for (k in seq_len(n_min)) { # until we have to sample from both groups
      # sample one value from the first group
      seen_XT <- hdp_XT_sampler_help(1, group_indx = indx_1st, seen_XT = seen_XT, c0 = c0, c = c, P00 = P00)

      # sample one value from the second group
      seen_XT <- hdp_XT_sampler_help(1, group_indx = indx_2nd, seen_XT = seen_XT, c0 = c0, c = c, P00 = P00)
    }

    # sample the remaining values from the group corresponding to the maximal number
    seen_XT <- hdp_XT_sampler_help(n_max - n_min, group_indx = max_indx, seen_XT = seen_XT, c0 = c0, c = c, P00 = P00)
  }
  return(seen_XT)
}
