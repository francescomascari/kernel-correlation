update_q_probs <- function(
    q_vec_k,
    l_tot,
    l_k,
    c = 1,
    c0 = 1) {
  # ---------------------------------------------------------------------------
  # Compute the un-normalized probabilities for the distribution of one customer
  # across the tables for the same dish, given the other customers' allocations.
  # Recall: the first entry corresponds to the probability to create a new table
  # ---------------------------------------------------------------------------
  #
  # Arguments:
  #   q_vec_k     : numeric vector – disposition of customers per table
  #                                  for the specified dish across restaurants.
  #
  #   l_tot       : numeric – number of tables across restaurants.
  #
  #   l_k         : numeric vector – number of tables for the
  #                                  specified dish across restaurants.
  #
  #   c0, c       : numeric – concentration parameters of the hDP (default 1).
  #
  # Returns:
  #   numeric vector – un-normalized probabilities for each new assignment

  # number of tables
  len_q <- length(q_vec_k)

  # if no tables left for the specified dish
  if (l_k == 0) {
    # return a vector, giving full probability to a new table
    return(c(1, rep(0, len_q)))
  }

  # inizialize the vector of un-normalized probabilities
  un_probs <- rep(numeric(0), len_q + 1)

  # assign the un-normalized probability of a new table in log scale
  un_probs[1] <- log(c) - log(c0 + l_tot) + log(l_k)

  for (i in seq_len(len_q)) {
    # assign the un-normalized probability of an existing table in log scale
    un_probs[i + 1] <- log(q_vec_k[i])
  }

  # normalize the probabilites
  probs <- exp(un_probs - logSumExp(un_probs))

  return(probs)
}