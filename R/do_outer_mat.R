do_outer_mat <- function(
    vals1,
    vals2,
    kernel = "gaussian",
    par_k = list(sigma = 1, beta = 1, left_lim = -Inf, right_lim = +Inf)) {
  # ---------------------------------------------------------------------------
  # Compute the matrix with entries given by k(x,y)
  # for the kernel k specified in `kernel` with parameters in `par_k`
  # as x and y varies in `vals1` and `vals2`, respectively.
  # ---------------------------------------------------------------------------
  #
  # Arguments:
  #   vals1       : numeric vector – first values to build the matrix on.
  #
  #   vals2       : numeric vector – second values to build the matrix on.
  #
  #   kernel      : character – name of the kernel in use.
  #
  #   par_k       : list – parameters of `kernel`.
  #                  sigma         – parameter of the gaussian kernel.
  #                  beta          – parameter of the laplace kernel.
  #                  left_lim      – left limit of the reference interval
  #                                  of the set-wise kernel.
  #                  right_lim     – right limit of the reference interval
  #                                  of the set-wise kernel.
  #
  # Returns:
  #   numeric matrix – the matrix given by k(x,y)
  #                    for x in `vals1` and y in `vals2`.

  # if `kernel` is gaussian
  if (kernel == "gaussian") {

    # build the matrix with the corresponding kernel
    quad_mat <- exp(-outer(vals1, vals2, "-")^2 / (2 * par_k$sigma^2))
  }

  # if `kernel` is laplace
  if (kernel == "laplace") {

    # build the matrix with the corresponding kernel
    quad_mat <- exp(-abs(outer(vals1, vals2, "-")) / par_k$beta)
  }

  # if `kernel` is setwise
  if (kernel == "setwise") {

    # compute the indicator functions for `vals1` and `vals2`
    ind1 <- (par_k$left_lim <= vals1) * (vals1 <= par_k$right_lim)
    ind2 <- (par_k$left_lim <= vals2) * (vals2 <= par_k$right_lim)

    # build the matrix with the corresponding kernel
    quad_mat <- outer(ind1, ind2)
  }

  # if `kernel` is linear
  if (kernel == "linear") {

    # build the matrix with the corresponding kernel
    quad_mat <- outer(vals1, vals2)
  }

  return(quad_mat)
}
