#### Generate data ####
if (FALSE) {
  
  N_vec <- 2*4^c(2,3,4,5)
  N_len <- length(N_vec)
  max_N <- max(N_vec)
  
  M_vec <- 5^c(2,3,4,5)
  M_len <- length(M_vec)
  max_M <- max(M_vec)
  
  set.seed(134)
  seen <- hDP_XT2(n1=max_N,n2=max_M,P00=runif,smpl_method = "alt",start = 1,c0=1,c=1)
  
  kernel_vec <- c("gaussian","laplace","setwise","linear")
  par_k = list(sigma = 1, beta = 1, set_left_lim = 0,set_right_lim = 0.95)
  
  
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
  for (i in 1:N_len*M_len) {
    N <- val_mat$N[i]
    M <- val_mat$M[i]
    
    seen_now <- rbind(seen[1:N,],seen[(max_N+1):(max_N+M),])
    seen_now_mat <- hDP_XT2_to_mat(seen_now)[,c(1,3,4)]
    
    for (k in kernel_vec) {
      sum_vals <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","do_diag_vec","quad.form","quad.3form","hDP_var_help_int_tabs","hDP_cov_help_int_tabs","logSumExp","update_q_probs","gibbs_tabs","gibbs_smpl_corr")) %dopar%
      {gibbs_anal_corr(N_baseln=N_baseln,seen = seen_now_mat,baseln=runif,c0=1,c=1,kernel = k,par_k = par_k)}
      val_mat[i,k] <- sum_vals/reps
    }
  }
  stopImplicitCluster()
  
  save(val_mat,file = "simulations/convergence_rate/convergence_rate.RData")
} else {
  load("simulations/convergence_rate/convergence_rate.RData")
}
filter(val_mt)

val_mat <- dplyr::filter(val_mat, N != 16 | M != 25)

my_mat <- rbind(val_mat[,1:2],val_mat[,1:2],val_mat[,1:2],val_mat[,1:2])
my_mat$corr <- c(val_mat$gaussian,val_mat$laplace,val_mat$linear,val_mat$setwise)
my_mat$kernel <- as.factor(c(rep("Gaussian",nrow(my_mat)/4),rep("Laplace",nrow(my_mat)/4),rep("Linear",nrow(my_mat)/4),rep("Setwise",nrow(my_mat)/4)))
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
  
