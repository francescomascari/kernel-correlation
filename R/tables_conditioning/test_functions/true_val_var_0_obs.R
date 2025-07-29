true_val_var_0_obs <- function(
  l_tot,
  n1_tot,
  n2_tot,
  group,
  c0,
  c,
  int_diag_baseln,
  int_cross_baseln,
  int_baseln_0,
  k_00){
    int_diag_post <- (c0*int_diag_baseln + l_tot*k_00)/(c0 + l_tot)
    int_cross_post <- (c0^2*int_cross_baseln + 2*c0*l_tot*int_baseln_0 + l_tot^2*k_00)/(c0 + l_tot)^2
    int_post_0 <- (c0*int_baseln_0 + l_tot*k_00)/(c0 + l_tot)
    
    if(group == 1){
      n_tot <- n1_tot
    }
    else{
      n_tot <- n2_tot
    }
    return((c^2*(c0+l_tot+c+n_tot+1)*(int_diag_post - int_cross_post)/(c0+l_tot+1) + c*n_tot*(int_diag_post - 2*int_post_0 + k_00))/((c+n_tot+1)*(c+n_tot)^2))
}