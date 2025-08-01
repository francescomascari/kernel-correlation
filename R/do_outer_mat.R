do_outer_mat <- function(
    vals1,
    vals2,
    kernel = "gaussian",
    par_k = list(sigma = 1, beta = 1, set_left_lim = -Inf,set_right_lim = +Inf)){
      if(kernel == "gaussian"){
        quad_mat <- exp(-outer(vals1,vals2,"-")**2/(2*par_k$sigma^2))
      }
      if(kernel == "laplace"){
        quad_mat <- exp(-abs(outer(vals1,vals2,"-"))/par_k$beta)
      }
      if(kernel == "setwise"){
        quad_mat <- outer((par_k$set_left_lim <= vals1)*(vals1 <= par_k$set_right_lim),(par_k$set_left_lim <= vals2)*(vals2 <= par_k$set_right_lim))
      }
      if(kernel == "linear"){
        quad_mat <- outer(vals1,vals2)
      }
      return(quad_mat)
}
