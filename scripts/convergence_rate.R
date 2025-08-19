#### Generate data ####
if (TRUE) {
  
  N_vec <- 4^c(2,3,4,5)
  max_N <- max(N_vec)
  
  M_vec <- 5^c(2,3,4,5)
  max_M <- max(M_vec)
  
  set.seed(1234)
  seen <- hdp_XT_sampler(n1=max_N,n2=max_M,smpl_method = "alt",start = 1,c0=1,c=1,P00=runif)
  
  kernel_vec <- c("gaussian","laplace","setwise","linear")
  par_k = list(sigma = 1, beta = 1, left_lim = 0, right_lim = 0.95)
  
  
  val_mat <- expand.grid(N=N_vec,M=M_vec)
  val_mat$gaussian<- NA
  val_mat$laplace <- NA
  val_mat$setwise <- NA
  val_mat$linear <- NA
  
  
  require("doParallel")
  n_Cores <- detectCores()
  cluster <- makeCluster(n_Cores-1)
  registerDoParallel(cluster)
  
  
  N_baseln <-10000
  reps <- 10
  
  nrow_mat <- length(N_vec)*length(M_vec)
  for (i in 1:nrow_mat) {
    N <- val_mat$N[i]
    M <- val_mat$M[i]
    
    seen_now <- rbind(seen[1:N,],seen[(max_N+1):(max_N+M),])
    seen_now_mat <- hDP_XT2mat(seen_now)[,c("X", "Ncusts1", "Ncusts2")]
    
    for (k in kernel_vec) {
      sum_vals <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","quad.form","quad.3form","logSumExp","update_q_probs","gibbs_tabs","hdp_corr_anlys")) %dopar%
      {hdp_corr_anlys(N_baseln=N_baseln,seen = seen_now_mat,baseln=runif,c0=1,c=1,kernel = k,par_k = par_k)}
      val_mat[i,k] <- sum_vals/reps
    }
  }
  stopImplicitCluster()
  
  save(val_mat, file = "output/results/convergence_rate.RData")
} else {
  load("output/convergence_rate.RData")
}

lm_data <- data.frame(slope = rep(NA,4),p.val = rep(NA,4),r.sq = rep(NA,4))
rownames(lm_data) <- c("Gaussian","Laplace","Setwise","Linear")

for(k in kernel_vec <- c("Gaussian","Laplace","Setwise","Linear")){
  model_lm <- lm(log(corr) ~ log(NM),data=dplyr::filter(my_mat,kernel==k))
  lm_data[k,] <- c(summary(model_lm)$coefficients["log(NM)",c(1,4)],summary(model_lm)$r.squared)
}

save(lm_data,file = "simulations/convergence_rate/convergence_rate_lm_data.RData")

my_mat <- rbind(val_mat[,1:2],val_mat[,1:2],val_mat[,1:2],val_mat[,1:2])
my_mat$corr <- c(val_mat$gaussian,val_mat$laplace,val_mat$linear,val_mat$setwise)
my_mat$kernel <- as.factor(c(rep("Gaussian",16),rep("Laplace",16),rep("Linear",16),rep("Setwise",16)))
my_mat$NM <- my_mat$N*my_mat$M

my_new_mat <- my_mat %>%
  group_by(kernel) %>%
  mutate(corr.smooth = exp(predict(lm(log(corr) ~ log(NM)))))

loadfonts()
ggplot(data = my_new_mat)+ 
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  geom_line(aes(x = N*M, y = corr.smooth,color = kernel),linewidth = 2) +
  geom_point(aes(x = N*M, y = corr,color = kernel, shape = kernel),size = 6) +
  #geom_point(aes(x = N*M, y = corr.smooth,color = kernel, shape = kernel),width = 2,size = 5) +
  labs(x = bquote(bold(n[1]*n[2])), y = "Correlation") +
  scale_color_manual(values = c("Gaussian"="purple3", "Laplace"="darkorange3","Linear"="lightgreen","Setwise"="cyan3")) + # associate each value of fill to a color
  scale_shape_manual(values = c("Gaussian"=15, "Laplace"=16,"Linear"=18,"Setwise"=17)) +
  theme_classic() +
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
        panel.grid = element_line(size = 1.5),
        legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(family="LM Roman 10",size=20,face="bold"))

  
