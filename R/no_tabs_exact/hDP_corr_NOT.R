hDP_corr_NOT <- function(
  baseln = runif,
  N_baseln = 10000,
  c = 1,
  c0 = 1,
  seen = matrix(numeric(0),ncol = 3),
  ...){
    baseln_smpl <- baseln(N_baseln)
    
    X_vec <- seen[,1]
    
    n1_vec <- seen[,2]
    n1_tot <- sum(n1_vec)
    
    n2_vec <- seen[,3]
    n2_tot <- sum(n2_vec)
    
    if(n1_tot != 0 || n2_tot != 0){
      l_df <- cond_law_tabs2_log(n1_vec = n1_vec,n2_vec = n2_vec,c0 = c0,c = c)
      
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
      mean_of_var1 <- mean_fun_ell(FUN = function(v){hDP_var_help_vec_tabs(int_cross_baseln,int_cross_custs1,int_cross_baseln_custs1,int_diag_baseln,int_diag_custs1,
                                                                           outer_mat_dishes,outer_mat_baseln_dishes,diag_vec_dishes,
                                                                           N_baseln,v,n1_vec,c = c,c0 = c0)},l_df=l_df)
      
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
      mean_of_var2 <- mean_fun_ell(FUN = function(v){hDP_var_help_vec_tabs(int_cross_baseln,int_cross_custs2,int_cross_baseln_custs2,int_diag_baseln,int_diag_custs2,
                                                                           outer_mat_dishes,outer_mat_baseln_dishes,diag_vec_dishes,
                                                                           N_baseln,v,n2_vec,c = c,c0 = c0)},l_df=l_df)
      
      mean_c0sq <- mean_fun_ell(FUN = function(v){(c0/(c0+sum(v)))^2},l_df = l_df)
      mean_c0 <- mean_fun_ell(FUN = function(v){c0/(c0+sum(v))},l_df = l_df)
      mean_c0l_int <- mean_fun_ell(FUN = function(v){c0/((c0+sum(v))^2*N_baseln)*sum(colSums(outer_mat_baseln_dishes)*v)},l_df = l_df)
      mean_l_int <- mean_fun_ell(FUN = function(v){1/((c0+sum(v))*N_baseln)*sum(colSums(outer_mat_baseln_dishes)*v)},l_df = l_df)
      mean_lsq_int <- mean_fun_ell(FUN = function(v){1/(c0+sum(v))^2*quad.form(outer_mat_dishes,v)},l_df = l_df)
      mean_l_vec <- mean_fun_ell(FUN = function(v){v/(c0+sum(v))},l_df = l_df)
      mean_l_vec_int <- quad.form(outer_mat_dishes,mean_l_vec)
    }
    else{
      return(hDP_corr_tabs(baseln = baseln,N_baseln = N_baseln,seen = matrix(numeric(0),ncol=4),c=c,c0=c0))
    }
    
    mean_of_cov <- mean_fun_ell(FUN = function(v){hDP_cov_help_vec_tabs(int_cross_baseln,int_diag_baseln,
                                                                        outer_mat_dishes,outer_mat_baseln_dishes,diag_vec_dishes,
                                                                        N_baseln,v,n1_vec,n2_vec,
                                                                        c = c,c0 = c0)},l_df=l_df)
    
    var_post_baseln <- (mean_c0sq - mean_c0^2)*int_cross_baseln + 2*(mean_c0l_int - mean_c0*mean_l_int) + mean_lsq_int - mean_l_vec_int
    
    var1 <- mean_of_var1 + c^2/((c+n1_tot)^2)*var_post_baseln
    var2 <- mean_of_var2 + c^2/((c+n2_tot)^2)*var_post_baseln
    
    cov <- mean_of_cov + c^2/((c+n1_tot)*(c+n2_tot))*var_post_baseln
    
    return(cov/sqrt(var1*var2))
}