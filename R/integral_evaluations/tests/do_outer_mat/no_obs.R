#### TEST 1: Accuracy of the integral for different set choices ####

# All real numbers
N_baseln <- 10000
baseln_smpl <- rnorm(N_baseln)
outer_mat_cross_baseln <- do_outer_mat(baseln_smpl,baseln_smpl,kernel = "setwise")

true_val1A <- 1
approx_val1A <- sum(outer_mat_cross_baseln)/N_baseln^2
print(c(true_val1A,approx_val1A))


# Positive real numbers
N_baseln <- 10000
baseln_smpl <- rnorm(N_baseln)
outer_mat_cross_baseln <- do_outer_mat(baseln_smpl,baseln_smpl,kernel = "setwise",par_k = list(set_left_lim = 0,set_right_lim = +Inf))

true_val1B <- 1/4
approx_val1B <- sum(outer_mat_cross_baseln)/N_baseln^2
print(c(true_val1B,approx_val1B))


# Central 95%
N_baseln <- 10000
baseln_smpl <- rnorm(N_baseln)
quant_val <- 0.95
lim <- qnorm((1 + quant_val)/2)
outer_mat_cross_baseln <- do_outer_mat(baseln_smpl,baseln_smpl,kernel = "setwise",par_k = list(set_left_lim = -lim,set_right_lim = lim))

true_val1C <- quant_val^2
approx_val1C <- sum(outer_mat_cross_baseln)/N_baseln^2
print(c(true_val1C,approx_val1C))


# Right-most 5%
N_baseln <- 10000
baseln_smpl <- rnorm(N_baseln)
quant_val <- 0.95
quant_val_1sided <- (1 + quant_val)/2
lim <- qnorm(quant_val_1sided)
outer_mat_cross_baseln <- do_outer_mat(baseln_smpl,baseln_smpl,kernel = "setwise",par_k = list(set_left_lim = lim,set_right_lim = +Inf))

true_val1D <- (1 - quant_val_1sided)^2
approx_val1D <- sum(outer_mat_cross_baseln)/N_baseln^2
c(true_val1D,approx_val1D)


#### TEST 2: Convergence of the integral as N increases #### TO DOOOOOO


if (FALSE) {
  N_vec <- c(10,30,100,300,1000,3000,10000)
  l_N <- length(N_vec)
  reps <- 1000
  
  val_mat2A <- matrix(numeric(0),nrow = l_N, ncol = reps)
  
  for (j in 1:reps) {
    last_N <- 0
    baseln_smpl <- numeric(0)
    for (i in 1:l_N) {
      N <- N_vec[i]
      delta_N <- N - last_N
      baseln_smpl <- c(baseln_smpl,rnorm(delta_N))
      outer_mat_cross_baseln <- do_outer_mat(baseln_smpl,baseln_smpl)
      
      val_mat2A[i,j] <- sum(outer_mat_cross_baseln)/N^2
      
      last_N <- N
    }
  }
  no_obs_T2A <- cbind(N=N_vec,True=true_val1A,t(apply(val_mat2A,MARGIN=1,FUN=summary)),Sd=apply(val_mat2A,MARGIN=1,FUN=sd))
  colnames(no_obs_T2A) <- c("N","True", "Min", "Q1", "Q2","Mean","Q3","Max","Sd")
  
  #write.table(no_obs_T2A, file = "integral_evaluations/tests/no_obs_T2A.csv", col.names = TRUE, row.names = TRUE)
} else {
  no_obs_T2A <- read.table("integral_evaluations/tests/do_outer_mat/no_obs_T2A.csv")
}

show(no_obs_T2A)

my_plot(no_obs_T2A)
