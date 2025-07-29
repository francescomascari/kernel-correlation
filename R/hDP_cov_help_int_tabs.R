hDP_cov_help_int_tabs <- function(
    int_cross_baseln,
    int_cross_tabs,
    int_cross_baseln_tabs,
    int_diag_baseln,
    int_diag_tabs,
    N_baseln,
    l_tot,
    n1_tot,
    n2_tot,
    c = 1,
    c0 = 1,
    ...){ 
      return(c^2/((c+n1_tot)*(c+n2_tot)*(c0+l_tot+1)*(c0+l_tot))*(c0*int_diag_baseln + l_tot*int_diag_tabs
             - (c0^2*int_cross_baseln + 2*c0*l_tot*int_cross_baseln_tabs + l_tot^2*int_cross_tabs)/(c0+l_tot)))
}