hDP_var_help_vec_tabs <- function(
    int_cross_baseln,
    int_cross_custs,
    int_cross_baseln_custs,
    int_diag_baseln,
    int_diag_custs,
    outer_mat_dishes,
    outer_mat_baseln_dishes,
    diag_vec_dishes,
    N_baseln,
    l_vec,
    n_vec,
    c = 1,
    c0 = 1,
    ...){
      l_tot <- sum(l_vec)
      n_tot <- sum(n_vec)
      
      if(l_tot != 0){
        int_cross_tabs <- quad.form(outer_mat_dishes,l_vec)/l_tot^2
        int_cross_baseln_tabs <- sum(colSums(outer_mat_baseln_dishes)*l_vec)/(N_baseln*l_tot)
        
        int_diag_tabs <- sum(diag_vec_dishes*l_vec)/l_tot
        
        if(n_tot != 0){
          int_cross_tabs_custs <- drop(quad.3form(outer_mat_dishes,l_vec,n_vec))/(l_tot*n_tot)
        }
        else{
          int_cross_tabs_custs <- 0
        }
      }
      else{
        int_cross_tabs <- 0
        int_cross_baseln_tabs <- 0
        int_cross_tabs_custs <- 0
        
        int_diag_tabs <- 0
      }
      
      return(hDP_var_help_int_tabs(int_cross_baseln,int_cross_custs,int_cross_tabs,int_cross_baseln_custs,int_cross_baseln_tabs,int_cross_tabs_custs,
                                   int_diag_baseln,int_diag_custs,int_diag_tabs,
                                   N_baseln,l_tot,n_tot,c = c,c0 = c0))
}