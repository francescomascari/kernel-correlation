hDP_var_tabs <- function(
    baseln = rnorm,
    N_baseln = 10000,
    seen = matrix(numeric(0),ncol = 4),
    c = 10,
    c0 = 10,
    group = 1,
    ...){
  baseln_smpl <- baseln(N_baseln)
  
  X_vec <- seen[,1]
  l_vec <- seen[,2]
  l_tot <- sum(l_vec)
  
  if(group == 1){
    n_vec <- seen[,3]
  }
  else if(group == 2){
    n_vec <- seen[,4]
  }
  n_tot <- sum(n_vec)
  
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
  
  return(hDP_var_help_vec_tabs(int_cross_baseln,int_cross_custs,int_cross_baseln_custs,int_diag_baseln,int_diag_custs,
                               outer_mat_dishes,outer_mat_baseln_dishes,diag_vec_dishes,
                               N_baseln,l_vec,n_vec,c = c,c0 = c0))
}