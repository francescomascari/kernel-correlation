do_diag_vec <- function(
    vals,
    kernel = "gaussian",
    par_k = list(sigma = 1, beta = 1, set_left_lim = -Inf,set_right_lim = +Inf)){
      if(kernel == "gaussian"){
        N <- length(vals)
        return(rep(1,N))
      }
      if(kernel == "laplace"){
        N <- length(vals)
        return(rep(1,N))
      }
      if (kernel == "setwise"){
        return((par_k$set_left_lim <= vals)*(vals <= par_k$set_right_lim))
      }
      if (kernel == "linear"){
        return(vals**2)
      }
}
