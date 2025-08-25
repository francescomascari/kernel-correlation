if(FALSE){
  set.seed(1234)
  seen <- hDP_XT2(n1=10,n2=10,P00=runif,smpl_method = "alt",start = 1,c0=1,c=1)
  seen_mat <- hDP_XT2_to_mat(seen)[,c(1,3,4),drop = FALSE]
} else {
  load("simulations/total_cov_VS_sampling/Poisson3_data.rda")
}

if(FALSE){
  require("doParallel")
  n_Cores <- detectCores()
  cluster <- makeCluster(n_Cores-1)
  registerDoParallel(cluster)
  
  reps <- 100
  
  smpl_vals <- foreach(i = seq_len(reps), .combine = c, .inorder = FALSE, .export = c("gibbs_smpl_corr","do_outer_diag","do_outer_mat","hDP_X2_mat_data","hDP_mat_X2_help","logSumExp","update_q_probs","gibbs_tabs")) %dopar%
    {gibbs_smpl_corr(baseln = runif,seen = seen_mat, M = 10000,c=1,c0=1)}
  anal_vals <- foreach(i = seq_len(reps), .combine = c, .inorder = FALSE, .export = c("gibbs_anal_corr","do_outer_diag","do_outer_mat","hDP_X2_mat_data","hDP_mat_X2_help","logSumExp","update_q_probs","gibbs_tabs")) %dopar%
    {gibbs_anal_corr(baseln = runif,seen = seen_mat,max_iter = 10,c=1,c0=1)}
  
  df_new <- data.frame(Value=c(smpl_vals,anal_vals),Method=as.factor(c(rep("Sampling",100),rep("Analytics",100))))
  save(df_new,file="smpl_vs_anal.RData")
  
} else {
  load("smpl_vs_anal.RData")
}

df_new <- data.frame(df_new)

quartiles_anal <- df_new %>%
  dplyr::filter(Method == "Analytics") %>% 
  reframe(y = round(quantile(Value, c(.25, .5, .75)),digits=4))

quartiles_smpl <- df_new %>% 
  dplyr::filter(Method == "Sampling") %>% 
  reframe(y = round(quantile(Value, c(.25, .5, .75)),digits=4))


options(digits = 6)
ggplot(df_new, aes(y = Value, x = Method)) +
  geom_boxplot(aes(fill = Method),linewidth=0.75,alpha = 0.4, show.legend = FALSE) +
  geom_text(data = quartiles_anal, aes(x = "Analytics", y = y + sign(y-y[2])*0.05, label = y),
            nudge_x = .4,
            hjust = 0,
            size = 5,
            family = "LM Roman 10") +
  geom_text(data = quartiles_smpl, aes(x = "Sampling", y = y + sign(y-y[2])*0.005, label = y),
            nudge_x = .4,
            hjust = 0,
            size = 5,
            family = "LM Roman 10") +
  scale_fill_manual(name = "Method", values = c("Sampling"="red3","Analytics"="blue3")) +
  labs(y = "Correlation", x = "Method") +
  theme_classic() +
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"))

save
