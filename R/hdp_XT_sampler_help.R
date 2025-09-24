hdp_XT_sampler_help <- function(
    n,
    group_indx,
    seen_XT = matrix(numeric(0),
                     ncol = 3,
                     dimnames = list(NULL, c("X", "T", "group"))),
    c0 = 1,
    c = 1,
    P00 = runif) {
  # ---------------------------------------------------------------------------
  # Draw n observations for a single group from the hierarchical DP (hDP) model.
  # If `seen_XT` is supplied, draws are conditional on the values in that matrix.
  # Sampling follows the Restaurant franchise scheme with table (cluster) indicators:
  #   * A new customer joins an existing table with prob. Ncusts / (c + Ncusts)
  #   * Otherwise, a new table is created.
  #       – The table serves an existing dish with prob. Ntabs / (c0 + Ntabs)
  #       – Otherwise, it introduces a brand-new dish from the baseline measure P00.
  # The function updates the shared `seen_XT` matrix.
  # ---------------------------------------------------------------------------
  #
  # Arguments
  #   n           : integer – number of customers to add to the given group.
  #
  #   group_indx  : integer (1 or 2) – which group these customers belong to.
  #
  #   seen_XT     : numeric matrix with columns,
  #                 (using the restaurant franchise metaphor)
  #                   X      – label of the dish
  #                   T      – label of the table
  #                   group  – group of the observation (1/2)
  #                 If empty (default), sampling is unconditional.
  #
  #   c, c0       : numeric – concentration parameters of the hDP.
  #
  #   P00         : function – baseline measure generator (default `runif`).
  #
  # Returns
  #   numeric matrix with columns : updated `seen_XT` matrix
  #                                 (same column structure as above).

  if (n == 0) { # if no value is requested to be sampled
    # return the original `seen_XT` matrix
    return(seen_XT)
  }

  # for the data in the `seen_XT` matrix, record
  row_seen <- nrow(seen_XT) # number of rows

  # create a matrix without duplicate tables
  # unique_seen_XT <- matrix(seen_XT[!duplicated(seen_XT[,"T"]), c("X","T")], ncol = 2, dimnames = list(NULL, c("X", "T")))
  unique_seen_XT <- seen_XT[!duplicated(seen_XT[, "T"]), c("X", "T"), drop = FALSE]
  n_tabs <- length(unique_seen_XT) # number of tables
  n_custs <- sum(seen_XT[, "group"] == group_indx) # number of costumers

  for (k in seq_len(n)) { # for every new value to be sampled
    if (runif(1) < c / (c + n_custs)) {
      # new table
      tab <- runif(1)
      if (runif(1) < c0 / (c0 + n_tabs)){
        # new dish
        dish <- P00(1)
      } else {
        # old dish
        dish <- sample(unique_seen_XT[, "X"], size = 1)
      }
      n_tabs <- n_tabs + 1

      seen_XT <- rbind(seen_XT, c(dish, tab, group_indx))
      unique_seen_XT <- rbind(unique_seen_XT, c(dish, tab))
    } else {
      # old table, old dish
      indx <- sample(seq_len(row_seen + k - 1), size = 1, prob = as.integer(seen_XT[,"group"] == group_indx))
      seen_XT <- rbind(seen_XT, seen_XT[indx, ])
    }
    n_custs <- n_custs + 1
  }
  return(seen_XT)
}
