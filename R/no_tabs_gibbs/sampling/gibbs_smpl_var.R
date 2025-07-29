gibbs_smpl_var <- function(
  baseln = runif,
  M = 10000,
  c = 1,
  c0 = 1,
  group = 1,
  seen = matrix(numeric(0),ncol = 3,dimnames = list(NULL,c("X", "Ncusts1", "Ncusts2"))),
  ...){
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
        if(group == 1){
          joint_vals <- hDP_X2_mat_data(2,0,seen = seen_mat,smpl_method = "alt",start = 1,c=c,c0=c0,P00=baseln)$new_X
          X1_vals <- c(X11_vals,joint_vals[joint_vals[,"group"]==1,1][1])
          X2_vals <- c(X12_vals,joint_vals[joint_vals[,"group"]==1,1][2])
        } else {
          joint_vals <- hDP_X2_mat_data(0,2,seen = seen_mat,smpl_method = "alt",start = 1,c=c,c0=c0,P00=baseln)$new_X
          X1_vals <- c(X11_vals,joint_vals[joint_vals[,"group"]==2,1][1])
          X2_vals <- c(X12_vals,joint_vals[joint_vals[,"group"]==2,1][2])
        }
      }
    } else {
      for(i in seq_len(M)){
        l_vec <- sapply(seen_q1_mat,length) + sapply(seen_q2_mat,length)
        seen_mat <- cbind(seen,"Ntabs"=l_vec) # we impose this table disposition
        
        if(group == 1){
          joint_vals <- hDP_X2_mat_data(2,0,seen = seen_mat,smpl_method = "alt",start = 1,c=c,c0=c0,P00=baseln)$new_X
          X1_vals <- c(X11_vals,joint_vals[joint_vals[,"group"]==1,1][1])
          X2_vals <- c(X12_vals,joint_vals[joint_vals[,"group"]==1,1][2])
        } else {
          joint_vals <- hDP_X2_mat_data(0,2,seen = seen_mat,smpl_method = "alt",start = 1,c=c,c0=c0,P00=baseln)$new_X
          X1_vals <- c(X11_vals,joint_vals[joint_vals[,"group"]==2,1][1])
          X2_vals <- c(X12_vals,joint_vals[joint_vals[,"group"]==2,1][2])
        }
  
        seen_q1_mat <- gibbs_tabs(l_vec,seen_q1_mat,c0=c0,c=c)
        l_vec <- sapply(seen_q1_mat,length) + sapply(seen_q2_mat,length)
        seen_q2_mat <- gibbs_tabs(l_vec,seen_q2_mat,c0=c0,c=c)
      }
    }
    
    return((sum(do_outer_diag(X1_vals,X2_vals,...)) - sum(do_outer_mat(X1_vals,X2_vals,...))/M)/(M-1))
}
