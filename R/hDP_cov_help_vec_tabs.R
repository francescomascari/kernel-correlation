hDP_cov_help_vec_tabs <- function(
    int_cross_baseln,
    int_diag_baseln,
    outer_mat_dishes,
    outer_mat_baseln_dishes,
    diag_vec_dishes,
    N_baseln,
    l_vec,
    n1_vec,
    n2_vec,
    c = 1,
    c0 = 1,
    ...){
      l_tot <- sum(l_vec)
      n1_tot <- sum(n1_vec)
      n2_tot <- sum(n2_vec)
      
      if(l_tot != 0){
        int_cross_tabs <- quad.form(outer_mat_dishes,l_vec)/l_tot^2
        int_cross_baseln_tabs <- sum(colSums(outer_mat_baseln_dishes)*l_vec)/(N_baseln*l_tot)
        
        int_diag_tabs <- sum(diag_vec_dishes*l_vec)/l_tot
      }
      else{
        int_cross_tabs <- 0
        int_cross_baseln_tabs <- 0
        
        int_diag_tabs <- 0
      }
      
      return(hDP_cov_help_int_tabs(int_cross_baseln,int_cross_tabs,int_cross_baseln_tabs,
                                   int_diag_baseln,int_diag_tabs,
                                   N_baseln,l_tot,n1_tot,n2_tot,c = c,c0 = c0))
}