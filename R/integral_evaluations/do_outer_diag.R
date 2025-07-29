do_outer_diag <- function(
    vals1,
    vals2,
    kernel = "gaussian",
    par_k = list(sigma = 1, beta = 1, set_left_lim = -Inf,set_right_lim = +Inf)){
  if(kernel == "gaussian"){
    vec_val <- exp(-(vals1-vals2)**2/(2*par_k$sigma^2))
  }
  if(kernel == "laplace"){
    vec_val <- exp(-abs(vals1-vals2)/par_k$beta)
  }
  if(kernel == "setwise"){
    vec_val <- (par_k$set_left_lim <= vals1)*(vals1 <= par_k$set_right_lim)*(par_k$set_left_lim <= vals2)*(vals2 <= par_k$set_right_lim)
  }
  if(kernel == "linear"){
    vec_val <- vals1*vals2
  }
  return(vec_val)
}
