hdp_sampler_help <- function(
    n,
    group_indx,
    seen = matrix(numeric(0),ncol=4,dimnames = list(NULL,c("X", "Ncusts1", "Ncusts2", "Ntabs"))),
    c0 = 1,
    c = 1,
    P00 = runif) {
  # ---------------------------------------------------------------------------
  # Draw n observations for a single group from the hierarchical DP (hDP) model.
  # If `seen` is supplied, draws are conditional on the counts in that matrix.
  # Sampling follows the Restaurant franchise scheme with table (cluster) indicators:
  #   * A new customer joins an existing table with prob. Ncusts / (c + Ncusts)
  #   * Otherwise a new table is created.
  #       – The table serves an existing dish with prob. Ntabs / (c0 + Ntabs)
  #       – Otherwise it introduces a brand-new dish from the baseline measure P00.
  # The function updates the shared `seen` count matrix and returns the updated
  # matrix plus the n newly generated (value, group) rows.
  # ---------------------------------------------------------------------------
  #
  # Arguments
  #   n           : integer – number of customers to add to the given group.
  #
  #   group_indx  : integer (1 or 2) – which group these customers belong to.
  #
  #   seen        : numeric matrix with columns
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
  #     seen   = updated `seen` matrix (same column structure as above),
  #     new_X  = matrix of newly sampled rows with columns:S
  #                X     – sampled value
  #                group – group index (matches `group_indx`)
  #   )

  # initialize the matrix of newly sampled values
  new_X <-  matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("X", "group")))

  if (n == 0) { # if no value is requested to be sampled
    # return the original `seen` matrix and the empty `new_X` matrix
    return(list(seen=seen,new_X=new_X))
  }

  # for the data in the `seen` matrix, record
  n_dishes <- nrow(seen) #number of unique values already sampled
  n_tabs <- sum(seen[, "Ntabs"]) #number of tables
  n_custs <- sum(seen[, group_indx + 1]) #number of costumers

  for (k in seq_len(n)) { # for every new value to be sampled
    if (runif(1) < c / (c + n_custs)) {
      # new table
      if (runif(1) < c0 / (c0 + n_tabs)) {
        # new dish
        dish <- P00(1)
        n_dishes <- n_dishes + 1

        # add the entry related to the new dish to the `seen` matrix
        new_dish <- c(dish, 0, 0, 1)
        new_dish[group_indx + 1] <- 1
        seen <- rbind(seen, new_dish)
      } else {
        # old dish

        # sample an index in seen with respect to the distribution of tables
        indx <- sample(1:n_dishes, size = 1, prob = seen[, "Ntabs"])

        # update the entry related to the sampled dish in the `seen` matrix
        dish <- seen[indx, "X"]
        seen[indx, "Ntabs"] <- seen[indx, "Ntabs"] + 1
        seen[indx, group_indx + 1] <- seen[indx, group_indx + 1] + 1
      }
      n_tabs <- n_tabs + 1
    } else {
      # old table, old dish

      # sample an index in seen with respect to the distribution of customers in the group
      indx <- sample(1:n_dishes, size = 1, prob = seen[, group_indx + 1])

      # update the entry related to the sampled dish in the `seen` matrix
      dish <- seen[indx, "X"]
      seen[indx, group_indx + 1] <- seen[indx, group_indx + 1] + 1
    }
    n_custs <- n_custs + 1

    # update the `new_X` matrix with the sampled value and the corresponding group
    new_X <- rbind(new_X, c(dish, group_indx))
  }

  return(list(seen = seen, new_X = new_X))
}