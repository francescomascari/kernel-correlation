#### Generate data ####
if (FALSE) {
  
  NM <- 1000
  
  set.seed(1234)
  seen <- hDP_XT2(n1=NM,n2=NM,P00=runif,smpl_method = "alt",start = 1,c0=1,c=1)
  seen_mat <- hDP_XT2_to_mat(seen)[,c(1,3,4)]
  
  require("doParallel")
  n_Cores <- detectCores()
  cluster <- makeCluster(n_Cores-1)
  registerDoParallel(cluster)
  
  N_baseln <- 10000
  reps <- 10
  
  #### Set-wise kernel ####
  
  right_lim_vec <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
  l_right_lim <- length(right_lim_vec)
  
  set_val <- c()
  for(i in seq_len(l_right_lim)){
    par_k <- list(set_left_lim = 0,set_right_lim = right_lim_vec[i])
    sum_res <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","do_diag_vec","quad.form","quad.3form","hDP_var_help_int_tabs","hDP_cov_help_int_tabs","logSumExp","update_q_probs","gibbs_tabs","gibbs_anal_corr")) %dopar%
      {gibbs_anal_corr(N_baseln=N_baseln,seen = seen_mat,baseln=runif,c0=1,c=1,kernel = "setwise",par_k = par_k)}
    
    set_val[i] <- sum_res/reps
  }
  
  
  #### Gaussian kernel ####
  sigma_vec <- c(1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4)
  l_sigma <- length(sigma_vec)
  
  gauss_val <- c()
  for(i in seq_len(l_sigma)){
    par_k <- list(sigma = sigma_vec[i])
    sum_res <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","do_diag_vec","quad.form","quad.3form","hDP_var_help_int_tabs","hDP_cov_help_int_tabs","logSumExp","update_q_probs","gibbs_tabs","gibbs_smpl_corr")) %dopar%
      {gibbs_anal_corr(N_baseln=N_baseln,seen = seen_mat,baseln=runif,c0=1,c=1,kernel = "gaussian",par_k = par_k)}
    
    gauss_val[i] <- sum_res/reps
  }
  
  
  #### Laplace kernel ####
  
  beta_vec <- c(1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4)
  l_beta <- length(beta_vec)
  
  laplace_val <- c()
  
  for(i in seq_len(l_beta)){
    par_k <- list(beta = beta_vec[i])
    sum_res <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","do_diag_vec","quad.form","quad.3form","hDP_var_help_int_tabs","hDP_cov_help_int_tabs","logSumExp","update_q_probs","gibbs_tabs","gibbs_smpl_corr")) %dopar%
      {gibbs_anal_corr(N_baseln=N_baseln,seen = seen_mat,baseln=runif,c0=1,c=1,kernel = "laplace",par_k = par_k)}
    
    laplace_val[i] <- sum_res/reps
  }
  
  stopImplicitCluster()
  
  
  df <- t(apply(rbind(set_val,gauss_val,laplace_val),MARGIN=1,FUN=summary))
  colnames(df) <- c("Min", "Q1", "Q2","Mean","Q3","Max")
  rownames(df) <- c("Setwise","Gaussian","Laplace")
  
  save(df,file = "simulations/kernel_stability/kernel_stability.RData")
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
