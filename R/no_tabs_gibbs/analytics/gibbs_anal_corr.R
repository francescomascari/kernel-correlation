gibbs_anal_corr <- function(
  baseln = runif,
  max_iter = 1000,
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
      seen_q1_mat <- lapply(n1_vec, FUN = function(times){rep(1,times)})
      seen_q2_mat <- lapply(n2_vec, FUN = function(times){rep(1,times)})
      
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
      
      var1 <- 0
      var2 <- 0
      cov <- 0
      valc0sq <- 0
      valc0 <- 0
      valc0l <- 0
      vall <- 0
      vallsq <- 0
      vall_vec <- 0
      
      for(i in 1:max_iter){
        l_vec <- sapply(seen_q1_mat,length) + sapply(seen_q2_mat,length)
        l_tot <- sum(l_vec)
        
        int_cross_tabs <- quad.form(outer_mat_dishes,l_vec)/l_tot^2
        int_cross_baseln_tabs <- sum(colSums(outer_mat_baseln_dishes)*l_vec)/(N_baseln*l_tot)
        int_diag_tabs <- sum(diag_vec_dishes*l_vec)/l_tot
        
        if(n1_tot!=0){
          int_cross_tabs_custs1 <- drop(quad.3form(outer_mat_dishes,l_vec,n1_vec))/(l_tot*n1_tot)
        }
        else{
          int_cross_tabs_custs1 <- 0
        }
        
        var1 <- var1 + hDP_var_help_int_tabs(int_cross_baseln,int_cross_custs1,int_cross_tabs,int_cross_baseln_custs1,int_cross_baseln_tabs,int_cross_tabs_custs1,
                                             int_diag_baseln,int_diag_custs1,int_diag_tabs,
                                             N_baseln,l_tot,n1_tot,c = c,c0 = c0)
        
        if(n2_tot!=0){
          int_cross_tabs_custs2 <- drop(quad.3form(outer_mat_dishes,l_vec,n2_vec))/(l_tot*n2_tot)
        }
        else{
          int_cross_tabs_custs2 <- 0
        }
        
        var2 <- var2 + hDP_var_help_int_tabs(int_cross_baseln,int_cross_custs2,int_cross_tabs,int_cross_baseln_custs2,int_cross_baseln_tabs,int_cross_tabs_custs2,
                                             int_diag_baseln,int_diag_custs2,int_diag_tabs,
                                             N_baseln,l_tot,n2_tot,c = c,c0 = c0)
        
        cov <- cov + hDP_cov_help_int_tabs(int_cross_baseln,int_cross_tabs,int_cross_baseln_tabs,
                                           int_diag_baseln,int_diag_tabs,
                                           N_baseln,l_tot,n1_tot,n2_tot,c = c,c0 = c0)
        
        valc0sq <- valc0sq + (c0/(c0+l_tot))^2
        valc0 <- valc0 + c0/(c0+l_tot)
        valc0l <- valc0l + c0*l_tot/(c0+l_tot)^2*int_cross_baseln_tabs
        vall <- vall + l_tot/(c0+l_tot)*int_cross_baseln_tabs
        vallsq <- vallsq + l_tot^2/(c0+l_tot)^2*int_cross_tabs
        vall_vec <- vall_vec + l_vec/(c0+l_tot)
        
        seen_q1_mat <- gibbs_tabs(l_vec,seen_q1_mat,c0=c0,c=c)
        l_vec <- sapply(seen_q1_mat,length) + sapply(seen_q2_mat,length)
        seen_q2_mat <- gibbs_tabs(l_vec,seen_q2_mat,c0=c0,c=c)
      }
    }
    else{
      return(hDP_corr_tabs(baseln = baseln,N_baseln = N_baseln,seen = matrix(numeric(0),ncol=4),c=c,c0=c0))
    }
    
    mean_valc0sq <- valc0sq/max_iter
    mean_valc0 <- valc0/max_iter
    mean_valc0l <- valc0l/max_iter
    mean_vall <- vall/max_iter
    mean_vallsq <- vallsq/max_iter
    mean_vall_vec <- vall_vec/max_iter
    
    var1_NOT <- var1/max_iter + c^2/(c+n1_tot)^2*((mean_valc0sq - (mean_valc0)^2)*int_cross_baseln + 2*mean_valc0l- 2*mean_valc0*mean_vall + mean_vallsq - quad.form(outer_mat_dishes,mean_vall_vec))
    var2_NOT <- var2/max_iter + c^2/(c+n2_tot)^2*((mean_valc0sq - (mean_valc0)^2)*int_cross_baseln + 2*mean_valc0l- 2*mean_valc0*mean_vall + mean_vallsq - quad.form(outer_mat_dishes,mean_vall_vec))
    cov_NOT <- cov/max_iter + c^2/((c+n1_tot)*(c+n2_tot))*((mean_valc0sq - (mean_valc0)^2)*int_cross_baseln + 2*mean_valc0l- 2*mean_valc0*mean_vall + mean_vallsq - quad.form(outer_mat_dishes,mean_vall_vec))
    
    return(cov_NOT/sqrt(var1_NOT*var2_NOT))
}
