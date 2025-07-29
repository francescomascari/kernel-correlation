true_val_cov_0_obs <- function(
  l_tot,
  n1_tot,
  n2_tot,
  c0,
  c,
  int_diag_baseln,
  int_cross_baseln,
  int_baseln_0,
  k_00){
    int_diag_post <- (c0*int_diag_baseln + l_tot*k_00)/(c0 + l_tot)
    int_cross_post <- (c0^2*int_cross_baseln + 2*c0*l_tot*int_baseln_0 + l_tot^2*k_00)/(c0 + l_tot)^2
    
    return(c^2*(int_diag_post - int_cross_post)/((c0+l_tot+1)*(c+n1_tot)*(c+n2_tot)))
}