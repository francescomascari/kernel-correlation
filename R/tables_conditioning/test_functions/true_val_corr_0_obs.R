true_val_corr_0_obs <- function(
  l_tot,
  n1_tot,
  n2_tot,
  ...){
    cov <- true_val_cov_0_obs(l_tot,n1_tot,n2_tot,...)
    var1 <- true_val_var_0_obs(l_tot,n1_tot,n2_tot,group=1,...)
    var2 <- true_val_var_0_obs(l_tot,n1_tot,n2_tot,group=2,...)
    return(cov/sqrt(var1*var2)) 
}