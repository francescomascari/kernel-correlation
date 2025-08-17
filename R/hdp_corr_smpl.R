hdp_corr_smpl <- function(
    baseln = runif,
    M = 10000,
    c = 1,
    c0 = 1,
    seen = matrix(numeric(0),
                  ncol = 3,
                  dimnames = list(NULL,c("X", "Ncusts1", "Ncusts2"))),
    ...) {
    
    X_vec <- seen[,1]
  
    n1_vec <- seen[,2]
    n1_tot <- sum(n1_vec)
    
    n2_vec <- seen[,3]
    n2_tot <- sum(n2_vec)
    
    if(n1_tot != 0 || n2_tot != 0){
      seen_q1_mat <- lapply(n1_vec, FUN = function(times){rep(1,times)})
      seen_q2_mat <- lapply(n2_vec, FUN = function(times){rep(1,times)})
    }
    
    X11_vals <- c()
    X12_vals <- c()
    X21_vals <- c()
    X22_vals <- c()
    if(n1_tot == 0 && n2_tot == 0){
      seen_mat <- matrix(numeric(0),ncol=4,dimnames = list(NULL,c("X", "Ncusts1", "Ncusts2", "Ntabs")))
      for(i in seq_len(M)){
        joint_vals <- hDP_X2_mat_data(2,2,seen = seen_mat,smpl_method = "alt",start = 1,c=c,c0=c0,P00=baseln)$new_X
        X11_vals <- c(X11_vals,joint_vals[joint_vals[,"group"]==1,1][1])
        X12_vals <- c(X12_vals,joint_vals[joint_vals[,"group"]==1,1][2])
        X21_vals <- c(X21_vals,joint_vals[joint_vals[,"group"]==2,1][1])
        X22_vals <- c(X22_vals,joint_vals[joint_vals[,"group"]==2,1][2])
      }
    } else {
      for(i in seq_len(M)){
        l_vec <- sapply(seen_q1_mat,length) + sapply(seen_q2_mat,length)
        seen_mat <- cbind(seen,"Ntabs"=l_vec) # we impose this table disposition
        
        joint_vals <- hDP_X2_mat_data(2,2,seen = seen_mat,smpl_method = "alt",start = 1,c=c,c0=c0,P00=baseln)$new_X
        X11_vals <- c(X11_vals,joint_vals[joint_vals[,"group"]==1,1][1])
        X12_vals <- c(X12_vals,joint_vals[joint_vals[,"group"]==1,1][2])
        X21_vals <- c(X21_vals,joint_vals[joint_vals[,"group"]==2,1][1])
        X22_vals <- c(X22_vals,joint_vals[joint_vals[,"group"]==2,1][2])
  
        seen_q1_mat <- gibbs_tabs(l_vec,seen_q1_mat,c0=c0,c=c)
        l_vec <- sapply(seen_q1_mat,length) + sapply(seen_q2_mat,length)
        seen_q2_mat <- gibbs_tabs(l_vec,seen_q2_mat,c0=c0,c=c)
      }
    }
    
    corr <- sum(do_outer_diag(X11_vals,X21_vals,...)) - sum(do_outer_mat(X11_vals,X21_vals,...))/M
    var1 <- sum(do_outer_diag(X11_vals,X12_vals,...)) - sum(do_outer_mat(X11_vals,X12_vals,...))/M
    var2 <- sum(do_outer_diag(X21_vals,X22_vals,...)) - sum(do_outer_mat(X21_vals,X22_vals,...))/M
    
    if(var1 >=0 && var2 >= 0){
      return(corr/sqrt(var1*var2))
    } else{
      return(NA)
    }
}
