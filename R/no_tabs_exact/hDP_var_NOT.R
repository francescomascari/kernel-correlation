hDP_var_NOT <- function(
  baseln = rnorm,
  N_baseln = 10000,
  c = 10,
  c0 = 10,
  seen = matrix(numeric(0),ncol = 3),
  group = 1,
  ...){
    baseln_smpl <- baseln(N_baseln)
    
    X_vec <- seen[,1]
    
    n1_vec <- seen[,2]
    n1_tot <- sum(n1_vec)
    
    n2_vec <- seen[,3]
    n2_tot <- sum(n2_vec)
    
    if(group == 1){
      n_vec <- n1_vec
      n_tot <- n1_tot
    }
    else{
      n_vec <- n2_vec
      n_tot <- n2_tot
    }
    if(n1_tot != 0 || n2_tot != 0){
      l_df <- cond_law_tabs2_log(n1_vec = n1_vec,n2_vec = n2_vec,c0 = c0,c = c)
      
      outer_mat_baseln <- do_outer_mat(vals1 = baseln_smpl,vals2 = baseln_smpl,...)
      outer_mat_dishes <- do_outer_mat(vals1 = X_vec,vals2 = X_vec,...)
      outer_mat_baseln_dishes <- do_outer_mat(vals1 = baseln_smpl,vals2 = X_vec,...)
      
      diag_vec_baseln <- do_diag_vec(vals = baseln_smpl,...)
      diag_vec_dishes <- do_diag_vec(vals = X_vec,...)
      
      int_cross_baseln <- sum(outer_mat_baseln)/N_baseln^2
      
      int_diag_baseln <- sum(diag_vec_baseln)/N_baseln
      
      if(n_tot != 0){
        int_cross_custs <- quad.form(outer_mat_dishes,n_vec)/n_tot^2
        int_cross_baseln_custs <- sum(colSums(outer_mat_baseln_dishes)*n_vec)/(N_baseln*n_tot)
        
        int_diag_custs <- sum(diag_vec_dishes*n_vec)/n_tot
      }
      else{
        int_cross_custs <- 0
        int_cross_baseln_custs <- 0
        
        int_diag_custs <- 0
      }
      
      mean_c0sq <- mean_fun_ell(FUN = function(v){(c0/(c0+sum(v)))^2},l_df = l_df)
      mean_c0 <- mean_fun_ell(FUN = function(v){c0/(c0+sum(v))},l_df = l_df)
      mean_c0l_int <- mean_fun_ell(FUN = function(v){c0/((c0+sum(v))^2*N_baseln)*sum(colSums(outer_mat_baseln_dishes)*v)},l_df = l_df)
      mean_l_int <- mean_fun_ell(FUN = function(v){1/((c0+sum(v))*N_baseln)*sum(colSums(outer_mat_baseln_dishes)*v)},l_df = l_df)
      mean_lsq_int <- mean_fun_ell(FUN = function(v){1/(c0+sum(v))^2*quad.form(outer_mat_dishes,v)},l_df = l_df)
      mean_l_vec <- mean_fun_ell(FUN = function(v){v/(c0+sum(v))},l_df = l_df)
      mean_l_vec_int <- quad.form(outer_mat_dishes,mean_l_vec)
    }
    else{
      return(hDP_var_tabs(baseln = baseln,N_baseln = N_baseln,seen = matrix(numeric(0),ncol=4),group=group))
    }
    mean_of_var <- mean_fun_ell(FUN = function(v){hDP_var_help_vec_tabs(int_cross_baseln,int_cross_custs,int_cross_baseln_custs,int_diag_baseln,int_diag_custs,
                                                                        outer_mat_dishes,outer_mat_baseln_dishes,diag_vec_dishes,
                                                                        N_baseln,v,n_vec,c = c,c0 = c0)},l_df=l_df)
    
    var_of_mean <- c^2/(c+n_tot)^2*((mean_c0sq - mean_c0^2)*int_cross_baseln + 2*(mean_c0l_int - mean_c0*mean_l_int) + mean_lsq_int - mean_l_vec_int)
    return(mean_of_var + var_of_mean)
}