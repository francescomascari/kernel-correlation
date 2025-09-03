hdp_mat_sampler_help <- function(
    n,
    group_indx,
    seen_mat = matrix(numeric(0),
                      ncol = 4,
                      dimnames = list(NULL, c("X", "Ncusts1", "Ncusts2", "Ntabs"))),
    c0 = 1,
    c = 1,
    P00 = runif) {
  # ---------------------------------------------------------------------------
  # Draw n observations for a single group from the hierarchical DP (hDP) model.
  # If `seen_mat` is supplied, draws are conditional on the counts in that matrix.
  # Sampling follows the Restaurant franchise scheme with table (cluster) indicators:
  #   * A new customer joins an existing table with prob. Ncusts / (c + Ncusts)
  #   * Otherwise a new table is created.
  #       – The table serves an existing dish with prob. Ntabs / (c0 + Ntabs)
  #       – Otherwise it introduces a brand-new dish from the baseline measure P00.
  # The function updates the shared `seen_mat` count matrix and returns the updated
  # matrix plus the n newly generated (value, group) rows.
  # ---------------------------------------------------------------------------
  #
  # Arguments
  #   n           : integer – number of customers to add to the given group.
  #
  #   group_indx  : integer (1 or 2) – which group these customers belong to.
  #
  #   seen_mat    : numeric matrix with columns
  #                 (using the restaurant franchise metaphor)
  #                   X        – unique dish values
  #                   Ncusts1  – frequency of customers in group 1
  #                   Ncusts2  – frequency of customers in group 2
  #                   Ntabs    – number of tables serving the dish
  #                 If empty, the process starts from scratch.
  #
  #   c, c0       : numeric – concentration parameters of the hDP.
  #
  #   P00         : function – baseline measure generator (default `runif`).
  #
  # Returns
  #   list(
  #     seen_mat   = updated `seen_mat` matrix (same column structure as above),
  #     new_X  = matrix of newly sampled rows with columns:S
  #                X     – sampled value
  #                group – group index (matches `group_indx`)
  #   )

  # initialize the matrix of newly sampled values
  new_X <-  matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("X", "group")))

  if (n == 0) { # if no value is requested to be sampled
    # return the original `seen_mat` matrix and the empty `new_X` matrix
    return(list(seen_mat = seen_mat, new_X = new_X))
  }

  # specify the group name depending on the group index
  group_name <- paste("Ncusts",group_indx,sep = "")

  # for the data in the `seen_mat` matrix, record
  n_dishes <- nrow(seen_mat) # number of unique values already sampled
  n_tabs <- sum(seen_mat[, "Ntabs"]) # number of tables
  n_custs <- sum(seen_mat[, group_name]) # number of costumers

  for (k in seq_len(n)) { # for every new value to be sampled
    if (runif(1) < c / (c + n_custs)) {
      # new table
      if (runif(1) < c0 / (c0 + n_tabs)) {
        # new dish
        dish <- P00(1)
        n_dishes <- n_dishes + 1

        # add the entry related to the new dish to the `seen_mat` matrix
        new_dish <- c(X = dish, Ncusts1 = 0, Ncusts2 = 0, Ntabs = 1)
        new_dish[group_name] <- 1
        seen_mat <- rbind(seen_mat, new_dish)
      } else {
        # old dish

        # sample an index in seen_mat with respect to the distribution of tables
        indx <- sample(1:n_dishes, size = 1, prob = seen_mat[, "Ntabs"])

        # update the entry related to the sampled dish in the `seen_mat` matrix
        dish <- seen_mat[indx, "X"]
        seen_mat[indx, "Ntabs"] <- seen_mat[indx, "Ntabs"] + 1
        seen_mat[indx, group_name] <- seen_mat[indx, group_name] + 1
      }
      n_tabs <- n_tabs + 1
    } else {
      # old table, old dish

      # sample an index in seen_mat with respect to the distribution of customers in the group
      indx <- sample(1:n_dishes, size = 1, prob = seen_mat[, group_name])

      # update the entry related to the sampled dish in the `seen_mat` matrix
      dish <- seen_mat[indx, "X"]
      seen_mat[indx, group_name] <- seen_mat[indx, group_name] + 1
    }
    n_custs <- n_custs + 1

    # update the `new_X` matrix with the sampled value and the corresponding group
    new_X <- rbind(new_X, c(dish, group_indx))
  }

  return(list(seen_mat = seen_mat, new_X = new_X))
}