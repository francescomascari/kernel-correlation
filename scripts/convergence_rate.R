## PART 1 : Data generation

# store all the different sample sizes for group 1
# and record the maximal
n1_all <- 4^c(2, 3, 4, 5)
n1_max <- max(n1_all)

# store all the different sample sizes for group 2
# and record the maximal
n2_all <- 5^c(2, 3, 4, 5)
n2_max <- max(n2_all)

# store the total number of cases
cases <- length(n1_all) * length(n2_all)

# set the seed for reproducibility
set.seed(1234)

# generate all the data up to the maximal sample sixes
# for both groups
seen <- hdp_XT_sampler(n1 = n1_max, n2 = n2_max, smpl_method = "alt", start = 1, c0 = 1, c = 1, P00 = runif)

# save the data matrix
write.csv(seen, file = "output/results/convergence_rate_seen.csv", row.names = FALSE)


## PART 2 : Computation of kernel correlation

# import the `seen` matrix from the previous step
seen <- as.matrix(read.csv("output/results/convergence_rate_seen.csv"))

# store the list of kernels under examination
# and their parameters
kernel_vec <- c("gaussian", "laplace", "setwise", "linear")
par_k <- list(sigma = 1, beta = 1, left_lim = 0, right_lim = 0.95)

# inizialize the matrix to store the values of the kernel correlation
val_mat <- expand.grid(n1 = n1_all, n2 = n2_all)
val_mat$gaussian <- NA
val_mat$laplace <- NA
val_mat$setwise <- NA
val_mat$linear <- NA

# set up parallelization
require("doParallel")
n_cores <- detectCores()
cluster <- makeCluster(n_Cores - 1)
registerDoParallel(cluster)

# set the number of repetitions to average the result on
reps <- 10

# for every case under study: each row of `val_mat``
for (i in seq_len(cases)) {
  # store the sample sizes under investigation
  n1 <- val_mat$n1[i]
  n2 <- val_mat$n2[i]

  # select the first `n1` observations from group 1
  # and the first `n2` observations from group 2
  seen_now <- rbind(seen[seen[, "group"] == 1, ][1:n1, ], seen[seen[, "group"] == 2, ][1:n2, ])

  # convert the observed data into the right format
  # for the `hdp_corr_anlys` function
  seen_now_mat <- hdp_XT2mat(seen_now)[, c("X", "Ncusts1", "Ncusts2"), drop = FALSE]

  # for every kernel under investigation
  for (k in kernel_vec) {
    # compute the kernel correlation `reps` times in parallel
    # and sum the realizations over the runs
    sum_vals <- foreach(i = seq_len(reps), .combine = sum, .inorder = FALSE, .export = c("do_outer_mat","quad.form","quad.3form","logSumExp","update_q_probs","gibbs_tabs","hdp_corr_anlys")) %dopar%
      {hdp_corr_anlys(seen = seen_now_mat, c0 = 1, c = 1, R = 1000, bsln = runif, M = 10000, kernel = k, par_k = par_k)}

    # store the mean value of the kernel correlation
    # across all the `reps` repetitions
    val_mat[i, k] <- sum_vals / reps
  }
}

# end parallelization
stopImplicitCluster()

# save the matrix with kernel correlation values
write.csv(val_mat, file = "output/results/convergence_rate_val_mat.csv", row.names = FALSE)


## PART 3 : Log-log plot for convergence rate

# import the `val_mat` matrix from the previous step as a data frame
val_mat <- read.csv("output/results/convergence_rate_val_mat.csv")

# pivot `val_mat` matrix to obtain a ggplot-friendly matrix
plot_mat <- pivot_longer(val_mat, cols = all_of(kernel_vec), names_to = "kernel", values_to = "corr") %>%
  # save the value of `n1 * n2` as a new column
  mutate(n1_n2 = n1 * n2) %>%
  # group it by `kernel`
  group_by(kernel) %>%
  # fit a linear model to infer the trend of the log correlation
  # with respect to the log of `n1 * n2`
  mutate(corr.smooth = exp(predict(lm(log(corr) ~ log(n1_n2)))))

# open the file to save the plot
png("output/plots/convergence_rate_plot.png", width = 900, height = 840, units = "px", res = 72)

# make a log log plot of the values of the kernel correlation vs `n1 * n2`
ggplot(data = plot_mat) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  # use the results of the linear model to plot the line
  geom_line(aes(x = n1_n2, y = corr.smooth, color = kernel), linewidth = 2) +
  # use the actual values to plot the points
  geom_point(aes(x = n1_n2, y = corr, color = kernel, shape = kernel), size = 6) +
  # add axis labels
  labs(x = bquote(bold(n[1] * n[2])), y = "Correlation") +
  # associate each kernel value to a fill color
  scale_color_manual(values = c("gaussian" = "purple3",
                                "laplace" = "darkorange3",
                                "linear" = "lightgreen",
                                "setwise" = "cyan3"),
                     labels = function(x) tools::toTitleCase(x)) +
  # associate each kernel value to a shape
  scale_shape_manual(values = c("gaussian" = 15,
                                "laplace" = 16,
                                "linear" = 18,
                                "setwise" = 17),
                     labels = function(x) tools::toTitleCase(x)) +
  # add a theme for better readability
  theme_classic() +
  # set the fonts of the plot to LaTeX style
  theme(axis.title = element_text(vjust = 0, family = "LM Roman 10", size = 30, face = "bold"),
        axis.text = element_text(size = 20, family = "LM Roman 10", face = "bold"),
        panel.grid = element_line(size = 1.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(family = "LM Roman 10", size = 20, face = "bold"))

# close the file to save the plot
dev.off()