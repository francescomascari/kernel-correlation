hDP_corr_tabs <- function(
  baseln = rnorm,
  N_baseln = 10000,
  seen = matrix(numeric(0),ncol = 4),
  c = 10,
  c0 = 10,
  ...){
    baseln_smpl <- baseln(N_baseln)
    
    X_vec <- seen[,1]
    l_vec <- seen[,2]
    l_tot <- sum(l_vec)
    
    n1_vec <- seen[,3]
    n1_tot <- sum(n1_vec)
    
    n2_vec <- seen[,4]
    n2_tot <- sum(n2_vec)
    
    outer_mat_baseln <- do_outer_mat(vals1 = baseln_smpl,vals2 = baseln_smpl,...)
    outer_mat_dishes <- do_outer_mat(vals1 = X_vec,vals2 = X_vec,...)
    outer_mat_baseln_dishes <- do_outer_mat(vals1 = baseln_smpl,vals2 = X_vec,...)
    
    diag_vec_baseln <- do_diag_vec(vals = baseln_smpl,...)
    diag_vec_dishes <- do_diag_vec(vals = X_vec,...)
    
    int_cross_baseln <- sum(outer_mat_baseln)/N_baseln^2
    
    int_diag_baseln <- sum(diag_vec_baseln)/N_baseln
    
    if(n1_tot != 0){
      int_cross_custs1 <- quad.form(outer_mat_dishes,n1_vec)/n1_tot^2
      int_cross_baseln_custs1 <- sum(colSums(outer_mat_baseln_dishes)*n1_vec)/(N_baseln*n1_tot)
      
      int_diag_custs1 <- sum(diag_vec_dishes*n1_vec)/n1_tot
    }
    else{
      int_cross_custs1 <- 0
      int_cross_baseln_custs1 <- 0
      
      int_diag_custs1 <- 0
    }
    
    if(n2_tot != 0){
      int_cross_custs2 <- quad.form(outer_mat_dishes,n2_vec)/n2_tot^2
      int_cross_baseln_custs2 <- sum(colSums(outer_mat_baseln_dishes)*n2_vec)/(N_baseln*n2_tot)
      
      int_diag_custs2 <- sum(diag_vec_dishes*n2_vec)/n2_tot
    }
    else{
      int_cross_custs2 <- 0
      int_cross_baseln_custs2 <- 0
      
      int_diag_custs2 <- 0
    }
    
    cov <- hDP_cov_help_vec_tabs(int_cross_baseln,int_diag_baseln,
                                 outer_mat_dishes,outer_mat_baseln_dishes,diag_vec_dishes,
                                 N_baseln,l_vec,n1_vec,n2_vec,c = c,c0 = c0)
    var1 <- hDP_var_help_vec_tabs(int_cross_baseln,int_cross_custs1,int_cross_baseln_custs1,int_diag_baseln,int_diag_custs1,
                                  outer_mat_dishes,outer_mat_baseln_dishes,diag_vec_dishes,
                                  N_baseln,l_vec,n1_vec,c = c,c0 = c0)
    var2 <- hDP_var_help_vec_tabs(int_cross_baseln,int_cross_custs2,int_cross_baseln_custs2,int_diag_baseln,int_diag_custs2,
                                  outer_mat_dishes,outer_mat_baseln_dishes,diag_vec_dishes,
                                  N_baseln,l_vec,n2_vec,c = c,c0 = c0)
    return(cov/sqrt(var1*var2))
}