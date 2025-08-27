## PART 1: Data generation

# set the seed for reproducibility
set.seed(1234)

# generate 10 data points for both groups
seen_mat <- hdp_mat_sampler(n1 = 10, n2 = 10, smpl_method = "alt", start = 1, c0 = 1, c = 1, P00 = runif)$seen_mat
seen <- seen_mat[, c("X", "Ncusts1", "Ncusts2"), drop = FALSE]

# save the data matrix
write.csv(seen, file = "output/results/smpl_vs_anlys_seen.csv", row.names = FALSE)


## PART 2 : Computation of kernel correlation

# import the `seen` matrix from the previous step
seen <- as.matrix(read.csv("output/results/smpl_vs_anlys_seen.csv"))

# set up parallelization
require("doParallel")
n_Cores <- detectCores()
cluster <- makeCluster(n_Cores - 1)
registerDoParallel(cluster)

# set the number of repetitions to build the box plots on
reps <- 100

# compute the kernel correlation `reps` times in parallel
# with the sampling-based algorithm and build a vector containing all the runs
smpl_vals <- foreach(i = seq_len(reps), .combine = c, .inorder = FALSE, .export = c("hdp_mat_sampler", "do_outer_mat", "logSumExp", "update_q_probs", "gibbs_tabs", "hdp_corr_smpl")) %dopar%
  {hdp_corr_smpl(seen = seen_mat, c0 = 1, c = 1, bsln = runif,  M = 10000)}

# compute the kernel correlation `reps` times in parallel
# with the analytics-based algorithm and build a vector containing all the runs
anal_vals <- foreach(i = seq_len(reps), .combine = c, .inorder = FALSE, .export = c("do_outer_mat", "quad.form", "quad.3form", "logSumExp", "update_q_probs", "gibbs_tabs", "hdp_corr_anlys")) %dopar%
  {hdp_corr_anlys(seen = seen_mat, c0 = 1, c = 1, R = 10, bsln = runif, M = 10000)}

# build a data frame storing all the values of the two different methods
df_val <- data.frame(Value = c(smpl_vals, anal_vals), Method = as.factor(c(rep("Sampling", 100), rep("Analytics", 100))))

# save the data matrix
write.csv(df_val, file = "output/results/smpl_vs_anlys_df_val.csv", row.names = FALSE)


## PART 3 : Plot of the boxplots for the repetitions of the two algorithm

# import the `df_val` matrix from the previous step
df_val <- read.csv("output/results/smpl_vs_anlys_df_val.csv")

# compute the quartiles for the values obtained for the sampling-based algorithm
quartiles_smpl <- df_val %>%
  dplyr::filter(Method == "Sampling") %>%
  reframe(y = round(quantile(Value, c(.25, .5, .75)), digits = 4))

# compute the quartiles for the values obtained for the analytics-based algorithm
quartiles_anlys <- df_val %>%
  dplyr::filter(Method == "Analytics") %>%
  reframe(y = round(quantile(Value, c(.25, .5, .75)), digits = 4))

# open the file to save the plot
png("output/plots/smpl_vs_anlys_plot.png", width = 900, height = 840, units = "px", res = 72)

# make a boxplot of the values obtained according to the two algorithms
# over the `reps` run
ggplot() +
  geom_boxplot(data = df_val, aes(x  = Method, y = Value, fill = Method), linewidth = 0.75, show.legend = FALSE) +
  # add labels referring to the values of the quartiles for each boxplot
  geom_text(data = quartiles_anlys, aes(x = "Analytics", y = y + sign(y-y[2])*0.05, label = y),
            nudge_x = .4,
            hjust = 0,
            size = 5,
            family = "LM Roman 10") +
  geom_text(data = quartiles_smpl, aes(x = "Sampling", y = y + sign(y-y[2])*0.005, label = y),
            nudge_x = .4,
            hjust = 0,
            size = 5,
            family = "LM Roman 10") +
  # set the color of the boxplot for each method
  scale_fill_manual(name = "Method", values = c("Analytics" = "firebrick", "Sampling" = "forestgreen")) +
  # se the axis labels
  labs(y = bquote(bold(Corr[k])), x = NULL) +
  # add a theme for better readability
  theme_classic() +
  # set the fonts of the plot to LaTeX style
  theme(axis.title = element_text(vjust = 0, family = "LM Roman 10", size = 30, face = "bold"),
        axis.text = element_text(size = 20, family = "LM Roman 10", face = "bold"))

# close the file to save the plot
dev.off()