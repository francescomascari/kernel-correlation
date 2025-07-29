hDP_var_help_int_tabs <- function(
  int_cross_baseln,
  int_cross_custs,
  int_cross_tabs,
  int_cross_baseln_custs,
  int_cross_baseln_tabs,
  int_cross_tabs_custs,
  int_diag_baseln,
  int_diag_custs,
  int_diag_tabs,
  N_baseln,
  l_tot,
  n_tot,
  c = 1,
  c0 = 1,
  ...){
    return((c*(c+n_tot)*(c0+l_tot+c+1)/((c0+l_tot+1)*(c0+l_tot))*(c0*int_diag_baseln+l_tot*int_diag_tabs)
      - c^2*(c0+l_tot+c+n_tot+1)/((c0+l_tot+1)*(c0+l_tot)^2)*(c0^2*int_cross_baseln+2*c0*l_tot*int_cross_baseln_tabs+l_tot^2*int_cross_tabs)
      - 2*c*n_tot/(c0+l_tot)*(c0*int_cross_baseln_custs+l_tot*int_cross_tabs_custs)
      + n_tot*(c+n_tot)*int_diag_custs
      - n_tot^2*int_cross_custs)/((c+n_tot+1)*(c+n_tot)^2))
}