#### Generate data ####
if (FALSE) {
  
  NM_vec <- c(0,10,100,1000)
  max_NM <- max(NM_vec)
  nrow_mat <- length(NM_vec)
  
  set.seed(1234)
  seen <- hDP_XT2(n1=max_NM,n2=max_NM,P00=runif,smpl_method = "alt",start = 1,c0=1,c=1)
  
  require("doParallel")
  n_Cores <- detectCores()
  cluster <- makeCluster(n_Cores-1)
  registerDoParallel(cluster)
  
  N_baseln <- 10000
  reps <- 10
  
  val_mat_stab <- data.frame(M = NM_vec,N = NM_vec)
  
  val_mat_stab$gaussian1 <- NA
  val_mat_stab$gaussian2 <- NA
  val_mat_stab$gaussian3 <- NA
  
  val_mat_stab$laplace1 <- NA
  val_mat_stab$laplace2 <- NA
  val_mat_stab$laplace3 <- NA
  
  val_mat_stab$setwise1 <- NA
  val_mat_stab$setwise2 <- NA
  val_mat_stab$setwise3 <- NA
  
  
  #### Setwise kernel ####
  right_lim_vec <- c(.1,.5,.9)
  l_right_lim <- length(right_lim_vec)
  
  for(j in seq_len(l_right_lim)){
    par_k <- list(set_left_lim = 0,set_right_lim = right_lim_vec[j])
    
    sum_res <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","do_diag_vec","quad.form","quad.3form","hDP_var_help_int_tabs","hDP_cov_help_int_tabs","logSumExp","update_q_probs","gibbs_tabs","gibbs_anal_corr")) %dopar%
      {gibbs_anal_corr(N_baseln=N_baseln,baseln=runif,c0=1,c=1,kernel = "setwise",par_k = par_k)}
    
    val_mat_stab[1,8+j] <- sum_res/reps
    
    for (i in 2:nrow_mat) {
      
      N <- val_mat_stab$N[i]
      M <- val_mat_stab$M[i]
      
      seen_now <- rbind(seen[seq_len(N),],seen[max_NM + seq_len(M),])
      seen_now_mat <- hDP_XT2_to_mat(seen_now)[,c(1,3,4),drop=FALSE]
      
      sum_res <- foreach(r = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","do_diag_vec","quad.form","quad.3form","hDP_var_help_int_tabs","hDP_cov_help_int_tabs","logSumExp","update_q_probs","gibbs_tabs","gibbs_anal_corr")) %dopar%
        {gibbs_anal_corr(N_baseln=N_baseln,seen = seen_now_mat,baseln=runif,c0=1,c=1,kernel = "setwise",par_k = par_k)}
      
      val_mat_stab[i,8+j] <- sum_res/reps
    }
  }
  
  #### Gaussian kernel ####
  sigma_vec <- c(1e-3,1e0,1e3)
  l_sigma <- length(sigma_vec)
  
  for(j in seq_len(l_sigma)){
    par_k <- list(sigma = sigma_vec[j])
    
    sum_res <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","do_diag_vec","quad.form","quad.3form","hDP_var_help_int_tabs","hDP_cov_help_int_tabs","logSumExp","update_q_probs","gibbs_tabs","gibbs_anal_corr")) %dopar%
      {gibbs_anal_corr(N_baseln=N_baseln,baseln=runif,c0=1,c=1,kernel = "gaussian",par_k = par_k)}
    
    val_mat_stab[1,2+j] <- sum_res/reps
    
    for (i in 2:nrow_mat) {
      
      N <- val_mat_stab$N[i]
      M <- val_mat_stab$M[i]
      
      seen_now <- rbind(seen[seq_len(N),],seen[max_NM + seq_len(M),])
      seen_now_mat <- hDP_XT2_to_mat(seen_now)[,c(1,3,4),drop=FALSE]
      
      sum_res <- foreach(r = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","do_diag_vec","quad.form","quad.3form","hDP_var_help_int_tabs","hDP_cov_help_int_tabs","logSumExp","update_q_probs","gibbs_tabs","gibbs_anal_corr")) %dopar%
        {gibbs_anal_corr(N_baseln=N_baseln,seen = seen_now_mat,baseln=runif,c0=1,c=1,kernel = "gaussian",par_k = par_k)}
      
      val_mat_stab[i,2+j] <- sum_res/reps
    }
  }
  
  #### Laplace kernel ####
  beta_vec <- c(1e-3,1e0,1e3)
  l_beta <- length(beta_vec)
  
  for(j in seq_len(l_beta)){
    par_k <- list(beta = beta_vec[j])
    
    sum_res <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","do_diag_vec","quad.form","quad.3form","hDP_var_help_int_tabs","hDP_cov_help_int_tabs","logSumExp","update_q_probs","gibbs_tabs","gibbs_anal_corr")) %dopar%
      {gibbs_anal_corr(N_baseln=N_baseln,baseln=runif,c0=1,c=1,kernel = "laplace",par_k = par_k)}
    
    val_mat_stab[1,5+j] <- sum_res/reps
    
    for (i in 2:nrow_mat) {
      
      N <- val_mat_stab$N[i]
      M <- val_mat_stab$M[i]
      
      seen_now <- rbind(seen[seq_len(N),],seen[max_NM + seq_len(M),])
      seen_now_mat <- hDP_XT2_to_mat(seen_now)[,c(1,3,4),drop=FALSE]
      
      sum_res <- foreach(r = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","do_diag_vec","quad.form","quad.3form","hDP_var_help_int_tabs","hDP_cov_help_int_tabs","logSumExp","update_q_probs","gibbs_tabs","gibbs_anal_corr")) %dopar%
        {gibbs_anal_corr(N_baseln=N_baseln,seen = seen_now_mat,baseln=runif,c0=1,c=1,kernel = "laplace",par_k = par_k)}
      
      val_mat_stab[i,5+j] <- sum_res/reps
    }
  }

  save(val_mat_stab,file = "simulations/kernel_stability_tab/kernel_stability_tab.RData")
} else {
  load("simulations/kernel_stability/poisson_data.rda")
}
df
min <- tibble(kernel = c("Setwise","Gaussian","Laplace"),y = round(df[,"Min"],digits=7))
max <- tibble(kernel = c("Setwise","Gaussian","Laplace"),y = round(df[,"Max"],digits=7))

options(digits=6)
col <- c("cyan3","purple3","darkorange3")
ggplot(data = df) +
  geom_errorbar(aes(x = rownames(df),y=Mean,ymin = Min,ymax = Max),linewidth=2, color=col,width = 0.25) +
  geom_text(data = min, aes(x = kernel, y = y, label = y),
            nudge_x = .15,
            hjust = 0,
            size = 5,
            family = "LM Roman 10") +
  geom_text(data = max, aes(x = kernel, y = y, label = y),
            nudge_x = .15,
            hjust = 0,
            size = 5,
            family = "LM Roman 10") +
  scale_y_continuous(trans='log10',) +
  labs(x = "", y = "Correlation") +
  theme_classic() +
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
        panel.grid = element_line(linewidth = 1.5),)
