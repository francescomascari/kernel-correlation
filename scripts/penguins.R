## PART 1 : Functions to fix parameters given the correlation
##          (These functions are in the Supplementary Material of the article)

# RECALL
# - `t_sq`   : marginal variance of the observables in the Gaussian example
# - `s_sq`   : variance of the likelihood in the Gaussian example
# - `tau_sq` : variance of the mean parameters in the Gaussian example
# - `rho`    : correlation of the mean parameters in the Gaussian example
# - `c0`,`c` : concentration parameters of the hierarchical Dirichlet Process (hDP)
# - `v`      : kernel variance a priori for each group in both cases
# - `xi`     : kernel correlation a priori in both cases
# - `sigma`  : parameter of the gaussian kernel

# value of `s_sq` as a function of `t_sq`, `v`, and `sigma`
compute_s_sq <- function(t_sq, v, sigma){
  return(sigma^2 * ((v + sqrt(sigma^2 / (2 * t_sq + sigma^2)))^-2 - 1) / 2)
}

# value of `tau_sq` as a function of `t_sq`, `v`, and `sigma`
compute_tau_sq <- function(t_sq, v, sigma){
  return(t_sq - sigma^2 * ((v + sqrt(sigma^2 / (2 * t_sq + sigma^2)))^-2 - 1) / 2)
}

# value of `rho` as a function of `xi`, `t_sq`, `v`, and `sigma`
compute_rho <- function(xi, t_sq, v, sigma){
  # numerator
  up <- t_sq - sigma^2 * ((v * xi + sqrt(sigma^2 / (2 * t_sq + sigma^2)))^-2 - 1) / 2

  # denominator
  down <- t_sq - sigma^2 * ((v + sqrt(sigma^2 / (2 * t_sq + sigma^2)))^-2 - 1) / 2

  return(up / down)
}

# value of `c0` as a function of `xi`, `t_sq`, `v`, and `sigma`
compute_c0 <- function(xi, t_sq, v, sigma){
  return((1 - sqrt(sigma^2 / (2 * t_sq + sigma^2))) / (v * xi) - 1)
}

# value of `c` as a function of `xi`, `t_sq`, `v`, and `sigma`
compute_c <- function(xi, t_sq, v, sigma){
  return(((1 - sqrt(sigma^2 / (2 * t_sq + sigma^2))) / v - 1) / (1 - xi))
}


## PART 2 : Data upload and kernel correlation computation

# set the values of `t_sq` and `v`
t_sq <- 2
v <- 1 / 4

# set the value of `sigma` accordingly
sigma <- sqrt(t_sq / (1 / (1 - v)^2 - 1))

# compute the value of `s_sq` and `tau_sq` accordingly
s_sq <- compute_s_sq(t_sq, v, sigma)
tau_sq <- compute_tau_sq(t_sq, v, sigma)

# Load the penguins dataset
library(palmerpenguins)
data(penguins)

# Clean the penguins dataset
penguins_clean <- penguins %>% drop_na()

# create unbalanced datsaset
# Set a seed for reproducibility, so you get the same 10 females every time
set.seed(40)

# Separate all the male penguins
males <- penguins_clean %>%
  dplyr::filter(sex == "male")

# Randomly sample exactly 10 female penguins
females_sample <- penguins_clean %>%
  dplyr::filter(sex == "female") %>%
  slice_sample(n = 5)

# Combine all males with the 10-female sample to create the unbalanced dataset
penguins_unbalanced <- bind_rows(males, females_sample)

# Check the new counts of males and females
penguins_unbalanced %>%
  count(sex)

# Calculate Ncusts1 (female occurrences) and Ncusts2 (male occurrences)
# from the new 'penguins_unbalanced' dataset
flipper_counts <- penguins_unbalanced %>%
  group_by(flipper_length_mm, sex) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = sex, values_from = count, values_fill = 0) %>%
  rename(X = flipper_length_mm, Ncusts1 = male, Ncusts2 = female)

# save vectors of counts for each group
n1_all <- flipper_counts$Ncusts1
n2_all <- flipper_counts$Ncusts2

# compute the total number of subjects for each group
n1 <- sum(n1_all)
n2 <- sum(n2_all)

# compute the means for each group
X1_bar <- flipper_counts %>%
  summarise(sum(X * Ncusts1) / sum(Ncusts1)) %>%
  as.numeric()
X2_bar <- flipper_counts %>%
  summarise(sum(X * Ncusts2) / sum(Ncusts2)) %>%
  as.numeric()

# store the size of the samples
M <- 10000

# store the values of the kernel correlation
xi_vec <- c(0.01, 0.50, 0.99)

# store the total number of cases
cases <- length(xi_vec)


## PART 2.1 : Gaussian example

# allocate the data frame to save the samples for the gaussian example
X_gau <- data.frame(X = numeric(0),
                    Group = factor(levels = c("Group 1", "Group 2")),
                    Corr = factor(levels = as.character(xi_vec)))

# for every value of the kernel correlation
for (i in seq_len(cases)) {
  # store the value of the kernel correlation for the case under investigation
  corr <- xi_vec[i]

  # compute the corresponding value of `rho`
  rho <- compute_rho(corr, t_sq, v, sigma)

  # compute the marginal means
  theta1 <- tau_sq * ((s_sq + n2 * tau_sq * (1 - rho^2)) * n1 * X1_bar + s_sq * rho * n2 * X2_bar) / (s_sq^2 + (n1 + n2) * s_sq * tau_sq + n1 * n2 * tau_sq^2 * (1 - rho^2))
  theta2 <- tau_sq * ((s_sq + n1 * tau_sq * (1 - rho^2)) * n2 * X2_bar + s_sq * rho * n1 * X1_bar) / (s_sq^2 + (n1 + n2) * s_sq * tau_sq + n1 * n2 * tau_sq^2 * (1 - rho^2))

  # compute the marginal variances
  Sigma11 <- s_sq + s_sq * tau_sq * (s_sq + n2 * tau_sq * (1 - rho^2)) / (s_sq^2 + (n1 + n2) * s_sq * tau_sq + n1 * n2 * tau_sq^2 * (1 - rho^2))
  Sigma22 <- s_sq + s_sq * tau_sq * (s_sq + n1 * tau_sq * (1 - rho^2)) / (s_sq^2 + (n1 + n2) * s_sq * tau_sq + n1 * n2 * tau_sq^2 * (1 - rho^2))

  # set the seed for reproducibility
  set.seed(2)

  # generate the samples
  X1_gau <- rnorm(M, mean = theta1, sd = sqrt(Sigma11))
  X2_gau <- rnorm(M, mean = theta2, sd = sqrt(Sigma22))

  # add the samples to the `X_gau` data frame
  X_gau <- rbind(X_gau,
                 data.frame(X = X1_gau,
                            Group = factor("Group 1",levels = c("Group 1", "Group 2")),
                            Corr = factor(as.character(corr),levels = as.character(xi_vec))),
                 data.frame(X = X2_gau,
                            Group = factor("Group 2",levels = c("Group 1", "Group 2")),
                            Corr = factor(as.character(corr),levels = as.character(xi_vec))))
}

# save the samples for the Gaussian example
write.csv(X_gau, file = "output/results/penguins_X_gau.csv", row.names = FALSE)


## PART 2.2 : Hierarchical Dirichlet Process
if (n1 != 0 || n2 != 0) {
  seen_q1_mat <- lapply(n1_all, FUN = function(times) {rep(1, times)})
  seen_q2_mat <- lapply(n2_all, FUN = function(times) {rep(1, times)})
}
# allocate the data frame to save the samples for the hDP
X_hdp <- data.frame(X = numeric(0),
                    Group = factor(levels = c("Group 1", "Group 2")),
                    Corr = factor(levels = as.character(xi_vec)))

# for every value of the kernel correlation
for (i in seq_len(cases)) {
  # store the value of the kernel correlation for the case under investigation
  corr <- xi_vec[i]

  # compute the corresponding value of `c0` and `c`
  c0 <- compute_c0(corr, t_sq, v, sigma)
  c <- compute_c(corr, t_sq, v, sigma)

  # set the seed for reproducibility
  set.seed(2)

  # allocate the vectors to store the samples
  X1_hdp <- c()
  X2_hdp <- c()

  # initialize by assigning one table per each customer for each group
  seen_q1_mat <- lapply(n1_all, FUN = function(times) {rep(1, times)})
  seen_q2_mat <- lapply(n2_all, FUN = function(times) {rep(1, times)})

  # if at least one group has some values
  if (n1 != 0 || n2 != 0) {
    # store the frequencies of unique values for the tables across both groups
    l_vec <- sapply(seen_q1_mat, length) + sapply(seen_q2_mat, length)
  } else {
    # assign an empty vector
    l_vec <- c()
  }

  for (i in seq_len(M)) {
    # impose this table disposition
    seen_mat <- cbind(flipper_counts, "Ntabs" = l_vec)

    # sample one value from the first group and save it in `X1_gau`
    joint_vals1 <- hdp_mat_sampler(1, 0, seen = seen_mat, smpl_method = "alt", start = 1, c0 = c0, c = c, P00 = function(n) {rnorm(n, mean = 200, sd = 10)})$new_X
    X1_hdp <- c(X1_hdp, joint_vals1[joint_vals1[, "group"] == 1, "X"][1])

    # sample one value from the second group and save it in `X2_gau`
    joint_vals2 <- hdp_mat_sampler(0, 1, seen = seen_mat, smpl_method = "alt", start = 2, c0 = c0, c = c, P00 = function(n) {rnorm(n, mean = 200, sd = 10)})$new_X
    X2_hdp <- c(X2_hdp, joint_vals2[joint_vals2[, "group"] == 2, "X"][1])

    # Gibbs update
    # if at least one group has some values
    if (n1 != 0 || n2 != 0) {
      # store the frequencies of unique values for the tables across both groups
      seen_q1_mat <- gibbs_tabs(l_vec, seen_q1_mat, c0 = c0, c = c)
      l_vec <- sapply(seen_q1_mat, length) + sapply(seen_q2_mat, length)
      seen_q2_mat <- gibbs_tabs(l_vec, seen_q2_mat, c0 = c0, c = c)
      l_vec <- sapply(seen_q1_mat, length) + sapply(seen_q2_mat, length)
    } else {
      # assign an empty vector
      l_vec <- c()
    }
  }

  X_hdp <- rbind(X_hdp,
                 data.frame(X = X1_hdp,
                            Group = factor("Group 1", levels = c("Group 1", "Group 2")),
                            Corr = factor(as.character(corr), levels = as.character(xi_vec))),
                 data.frame(X = X2_hdp,
                            Group = factor("Group 2", levels = c("Group 1", "Group 2")),
                            Corr = factor(as.character(corr), levels = as.character(xi_vec))))

}

# save the samples for the hDP
write.csv(X_hdp, file = "output/results/penguins_X_hdp.csv", row.names = FALSE)


## PART 5 : Plots of the mean posterior predictive distributions

# import the `X_gau` and `X_hdp` data frames from the previous step
X_gau <- read.csv("output/results/penguins_X_gau.csv")
X_hdp <- read.csv("output/results/penguins_X_hdp.csv")

# for every value of the kernel correlation
for (i in seq_len(cases)) {
  # store the value of the kernel correlation for the case under investigation
  corr <- xi_vec[i]

  # plot the distribution of the sample corresponding to the value of
  # kernel correlation under investigation for the Gaussian example
  p1 <- ggplot(subset(X_gau, Corr == as.character(corr)), aes(x = X, color = Group, linetype = Group)) +
    # make two density plots: one for the inner gray part, one for the border
    geom_density(position = "identity", fill = scales::alpha("darkgray", 0.67), linewidth = 0) +
    geom_density(position = "identity", fill = NA, linewidth = 3) +
    # set the color and the linetype according to the group
    scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "forestgreen")) +
    scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "dashed")) +
    # add a theme for better readability
    theme_classic() +
    # set the fonts of the plot to LaTeX style
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 30),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 30, face = "bold")) +
    # set the limits of the axes
    xlim(c(175, 210)) +
    ylim(c(0, 0.75))

  # save the plot as a png file
  fancy_png(plot = p1, out_path = paste("output/plots/penguins_gau", corr, ".png", sep = ""))


  # plot the distribution of the sample corresponding to the value of
  # kernel correlation under investigation for the hDP
  p2 <- ggplot(subset(X_hdp, Corr == as.character(corr)), aes(x = X, y = after_stat(density), color = Group, linetype = Group)) +
    # make two histogram plots: one for the inner gray part, one for the border
    geom_histogram(position = "identity", fill = scales::alpha("darkgray", 0.67), linewidth = 0) + 
    geom_histogram(position = "identity", fill = NA, linewidth = 3) +
    # set the color and the linetype according to the group
    scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "forestgreen")) +
    scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "dashed")) +
    # add a theme for better readability
    theme_classic() +
    # set the fonts of the plot to LaTeX style
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 30),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 30, face = "bold")) +
    # set the limits of the axes
    xlim(c(160, 240)) +
    ylim(c(0, 0.2))

    # save the plot as a png file
    fancy_png(plot = p2, out_path = paste("output/plots/penguins_hdp", corr, ".png", sep = ""))
}


## PART 6: Tile plot of the absolute distances of the means
##         between the mean posterior predictive distributions
##         of the second group

# allocate the data frame to store all couples of kernel correlation values
corr_df <- expand.grid(corr1 = as.factor(xi_vec), corr2 = as.factor(xi_vec))
cases_cross <- cases^2

# inizialize the vectors to store the absolute distances
val_gauVSgau <- numeric(cases_cross)
val_hdpVShdp <- numeric(cases_cross)
val_gauVShdp <- numeric(cases_cross)

# for every coupe of kernel correlation values
for (i in seq_len(cases_cross)){

  # extract the values of kernel correlations under investigation
  corr1 <- corr_df$corr1[i]
  corr2 <- corr_df$corr2[i]

  # compute the absolute distances of the means
  # and save them in the vectors
  val_gauVSgau[i] <- abs(mean(subset(X_gau, Corr == as.character(corr1) & Group == "Group 2")$X) - mean(subset(X_gau, Corr == as.character(corr2) & Group == "Group 2")$X))
  val_hdpVShdp[i] <- abs(mean(subset(X_hdp, Corr == as.character(corr1) & Group == "Group 2")$X) - mean(subset(X_hdp, Corr == as.character(corr2) & Group == "Group 2")$X))
  val_gauVShdp[i] <- abs(mean(subset(X_gau, Corr == as.character(corr1) & Group == "Group 2")$X) - mean(subset(X_hdp, Corr == as.character(corr2) & Group == "Group 2")$X))
}

# add the values of the kernel correlation to the right data frame
gauVSgau_df <- cbind(corr_df, value = val_gauVSgau)
hdpVShdp_df <- cbind(corr_df, value = val_hdpVShdp)
gauVShdp_df <- cbind(corr_df, value = val_gauVShdp)

# store the maximal value of the kernel correlation
# to set the upper limit of the color gradient
max_val <- max(c(val_gauVSgau, val_hdpVShdp, val_gauVShdp))

# make a tile plot comparing cases within the Gaussian case
p <- ggplot(gauVSgau_df, aes(x = corr1, y = corr2, fill = value)) +
  geom_tile(color = "white", linewidth = 2) +
  # add the actual values of the absolute distances
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 10, fontface = "bold") +
  # add a gradient scale for readability
  scale_fill_gradientn(colours = c("forestgreen", "firebrick"),
                       limits = c(0, max_val)) +
  # add axis labels
  labs(x = "$\\mathbb{C}\\mathrm{orr}_{k}$: Gaussian model",
       y = "$\\mathbb{C}\\mathrm{orr}_{k}$: Gaussian model") +
  # add a theme for better readability
  theme_classic() +
  # set the fonts of the plot to LaTeX style
  theme(
    axis.title = element_text(vjust = 0, size = 40),
    axis.text  = element_text(size = 30),
    legend.position = "none"
  )

# save the plot as a png file
fancy_png(plot = p, out_path = "output/plots/penguins_gauVSgau.png")


# make a tile plot comparing cases within the hDP
p <- ggplot(hdpVShdp_df, aes(x = corr1, y = corr2, fill = value)) +
  geom_tile(color = "white", linewidth = 2) +
  # add the actual values of the absolute distances
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 10, fontface = "bold") +
  # add a gradient scale for readability
  scale_fill_gradientn(colours = c("forestgreen", "firebrick"),
                       limits = c(0, max_val)) +
  # add axis labels
  labs(x = "$\\mathbb{C}\\mathrm{orr}_{k}$: hDP",
       y = "$\\mathbb{C}\\mathrm{orr}_{k}$: hDP") +
  # add a theme for better readability
  theme_classic() +
  # set the fonts of the plot to LaTeX style
  theme(
    axis.title = element_text(vjust = 0, size = 40),
    axis.text  = element_text(size = 30),
    legend.position = "none"
  )

# save the plot as a png file
fancy_png(plot = p, out_path = "output/plots/penguins_hdpVShdp.png")


# make a tile plot comparing cases within the hDP
p <- ggplot(gauVShdp_df, aes(x = corr1, y = corr2, fill = value)) +
  geom_tile(color = "white", linewidth = 2) +
  # add the actual values of the absolute distances
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 10, fontface = "bold") +
  # add a gradient scale for readability
  scale_fill_gradientn(colours = c("forestgreen", "firebrick"),
                       limits = c(0, max_val)) +
  # add axis labels
  labs(x = "$\\mathbb{C}\\mathrm{orr}_{k}$: Gaussian case",
       y = "$\\mathbb{C}\\mathrm{orr}_{k}$: hDP") +
  # add a theme for better readability
  theme_classic() +
  # set the fonts of the plot to LaTeX style
  theme(
    axis.title = element_text(vjust = 0, size = 40),
    axis.text  = element_text(size = 30),
    legend.position = "none"
  )

# save the plot as a png file
fancy_png(plot = p, out_path = "output/plots/penguins_gauVShdp.png")
