## PART 1: Data generation

# store all the different sample sizes for both groups
# and record the maximal
n_all <- c(0, 10, 100, 1000)
n_max <- max(n_all)

# store the total number of cases
cases <- length(n_all)

# set the seed for reproducibility
set.seed(1234)

# generate all the data up to the maximal sample sixes
# for both groups
seen <- hdp_XT_sampler(n1 = n_max, n2 = n_max, P00 = runif, smpl_method = "alt", start = 1, c0 = 1, c = 1)

# save the data matrix
write.csv(seen, file = "output/results/kernel_stability_seen.csv", row.names = FALSE)


## PART 2: Computation of kernel correlation

# inizialize the the matrix to store the values of the kernel correlation?
val_mat <- data.frame(n1 = n_all, n2 = n_all)
val_mat$gaussian1 <- NA
val_mat$gaussian2 <- NA
val_mat$gaussian3 <- NA
val_mat$laplace1 <- NA
val_mat$laplace2 <- NA
val_mat$laplace3 <- NA
val_mat$setwise1 <- NA
val_mat$setwise2 <- NA
val_mat$setwise3 <- NA

# set up parallelization
require("doParallel")
n_Cores <- detectCores()
cluster <- makeCluster(n_Cores - 1)
registerDoParallel(cluster)

# set the number of repetitions to average the result on
reps <- 10


## PART 2.1: Gaussian kernel

# store the kernel name
k <- "gaussian"

# inizialize the different values of `sigma` for the Gaussian kernel
sigma_vec <- c(1e-3, 1e0, 1e3)
sigma_cases <- length(sigma_vec)

# for every value of sigma
for (j in seq_len(sigma_cases)) {

  # set the parameter `sigma` of the kernel correlation
  par_k <- list(sigma = sigma_vec[j])
  for (i in seq_len(cases)) {
    # store the sample sizes under investigation
    n1 <- val_mat$n1[i]
    n2 <- val_mat$n2[i]

    # select the first `n1` observations from group 1
    # and the first `n2` observations from group 2
    seen_now <- rbind(seen[seen[, "group"] == 1, ][seq_len(n1), ], seen[seen[, "group"] == 2, ][seq_len(n2), ])

    # convert the observed data into the right format
    # for the `hdp_corr_anlys` function
    seen_now_mat <- hdp_XT2mat(seen_now)[, c("X", "Ncusts1", "Ncusts2"), drop = FALSE]

    # compute the kernel correlation `reps` times in parallel
    # and sum the realizations over the runs
    sum_vals <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","quad.form","quad.3form","logSumExp","update_q_probs","gibbs_tabs","hdp_corr_anlys")) %dopar%
      {hdp_corr_anlys(seen = seen_now_mat, c0 = 1, c = 1, R = 1000, bsln = runif, M = 10000, kernel = k, par_k = par_k)}
    
    # store the mean value of the kernel correlation
    # across all the `reps` repetitions
    val_mat[i, paste(k, j, sep = "")] <- sum_vals / reps
  }
}


## PART 2.2: Laplace kernel

# store the kernel name
k <- "laplace"

# inizialize the different values of `beta` for the Laplace kernel
beta_vec <- c(1e-3, 1e0, 1e3)
beta_cases <- length(beta_vec)

# for every value of beta
for (j in seq_len(beta_cases)) {
  
  # set the parameter `beta` of the kernel correlation
  par_k <- list(beta = beta_vec[j])
  for (i in seq_len(cases)) {
    # store the sample sizes under investigation
    n1 <- val_mat$n1[i]
    n2 <- val_mat$n2[i]
    
    # select the first `n1` observations from group 1
    # and the first `n2` observations from group 2
    seen_now <- rbind(seen[seen[, "group"] == 1, ][seq_len(n1), ], seen[seen[, "group"] == 2, ][seq_len(n2), ])
    
    # convert the observed data into the right format
    # for the `hdp_corr_anlys` function
    seen_now_mat <- hdp_XT2mat(seen_now)[, c("X", "Ncusts1", "Ncusts2"), drop = FALSE]
    
    # compute the kernel correlation `reps` times in parallel
    # and sum the realizations over the runs
    sum_vals <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","quad.form","quad.3form","logSumExp","update_q_probs","gibbs_tabs","hdp_corr_anlys")) %dopar%
      {hdp_corr_anlys(seen = seen_now_mat, c0 = 1, c = 1, R = 1000, bsln = runif, M = 10000, kernel = k, par_k = par_k)}
    
    # store the mean value of the kernel correlation
    # across all the `reps` repetitions
    val_mat[i, paste(k, j, sep = "")] <- sum_vals / reps
  }
}


## PART 2.3: Set-wise kernel

# store the kernel name
k <- "setwise"

# inizialize the different values of `right_lim` for the set-wise kernel
right_lim_vec <- c(.1,.5,.9)
right_lim_cases <- length(right_lim_vec)

# for every value of `right_lim`
for (j in seq_len(right_lim_cases)) {
  
  # set the parameter `beta` of the kernel correlation
  par_k <- list(left_lim = 0, right_lim = right_lim_vec[j])
  for (i in seq_len(cases)) {
    # store the sample sizes under investigation
    n1 <- val_mat$n1[i]
    n2 <- val_mat$n2[i]
    
    # select the first `n1` observations from group 1
    # and the first `n2` observations from group 2
    seen_now <- rbind(seen[seen[, "group"] == 1, ][seq_len(n1), ], seen[seen[, "group"] == 2, ][seq_len(n2), ])
    
    # convert the observed data into the right format
    # for the `hdp_corr_anlys` function
    seen_now_mat <- hdp_XT2mat(seen_now)[, c("X", "Ncusts1", "Ncusts2"), drop = FALSE]
    
    # compute the kernel correlation `reps` times in parallel
    # and sum the realizations over the runs
    sum_vals <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","quad.form","quad.3form","logSumExp","update_q_probs","gibbs_tabs","hdp_corr_anlys")) %dopar%
      {hdp_corr_anlys(seen = seen_now_mat, c0 = 1, c = 1, R = 1000, bsln = runif, M = 10000, kernel = k, par_k = par_k)}
    
    # store the mean value of the kernel correlation
    # across all the `reps` repetitions
    val_mat[i, paste(k, j, sep = "")] <- sum_vals / reps
  }
}

# end parallelization
stopImplicitCluster()

# save the matrix with kernel correlation values
write.csv(val_mat, file = "output/results/kernel_stability_val_mat.csv", row.names = FALSE)
