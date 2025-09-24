gibbs_tabs <- function(
    l_vec,
    seen_q_vec,
    c0 = 1,
    c = 1) {
  # ---------------------------------------------------------------------------
  # Update the disposition of the latent variables ("tables")
  # in `seen_q_vec` according to the augmented Gibbs sampling scheme
  # for the augmented hierarchical DP (hDP).
  # ---------------------------------------------------------------------------
  #
  # Arguments:
  #   l_vec       : numeric vector – number of tables per dish
  #                                  across restaurants.
  #
  #   seen_q_vec  : list of numeric vectors – disposition of customers per table
  #                                           across restaurants.
  #
  #   c0, c       : numeric – concentration parameters of the hDP (default 1).
  #
  # Returns:
  #   list of numeric vectors – updated `seen_q_vec`.

  # store the total number of dishes
  n_dishes <- length(l_vec)
  # store the number of tables per dish
  l_tot <- sum(l_vec)

  # for every dish
  for (n_indx in seq_len(n_dishes)) {
    # select the customers per table for the specified dish
    start_q_vec <- seen_q_vec[n_indx][[1]]
    q_vec <- start_q_vec

    # consider the number of tables for the specified dish
    len_q <- length(start_q_vec)

    # for every table for the specified dish
    for (q_indx in seq_len(len_q)) {
      # select the number of customers for the specified table
      q <- start_q_vec[q_indx]

      # for every customer for the specified table
      for (i in seq_len(q)) {
        # if there is only one customer in the specified table
        if (q_vec[q_indx] == 1) {

          # update `l_vec`, removing the table
          l_vec[n_indx] <- l_vec[n_indx] - 1
          l_tot <- l_tot - 1
        }

        # remove the customer from the specified table
        q_vec[q_indx] <- q_vec[q_indx] - 1

        # compute the un-normalized probability distribution of the Gibbs update
        # for the resulting disposition of customers and tables for the specified dish
        probs <- update_q_probs(q_vec, l_tot, l_vec[n_indx], c0 = c0, c = c)

        # select the index of an existing table to assign the customer to
        # according to the probabilities specified in `probs`
        # recall: 0 corresponds to a new table being created
        new_indx <- sample(0:length(q_vec), size = 1, prob = probs)

        # if we assign the customer to a new table
        if (new_indx == 0) {
          # add an entry in `q_vec`
          q_vec <- c(q_vec, 1)

          # update the count in l_vec and l_tot
          l_vec[n_indx] <- l_vec[n_indx] + 1
          l_tot <- l_tot + 1
        } else {
          # update the entry in `q_vec`
          q_vec[new_indx] <- q_vec[new_indx] + 1
        }
      }
    }
    # update `seen_q_vec`
    seen_q_vec[n_indx][[1]] <- q_vec[q_vec != 0]
  }
  return(seen_q_vec)
}
