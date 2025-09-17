## PART 1 : Data generation

# store the sample size
M <- 1000

# set the seed for reproducibility
set.seed(1234)

# set the correlation of the model to 1
# build the variance-covariance matrix
rho <- 1
var_cov_mat <- matrix(c(1, rho, rho, 1), nrow = 2)

# generate `n` independent realization of thetas
thetas <- rmvnorm(M, sigma = var_cov_mat)

# generate one independent sample for each group for each realization of thetas
X <- rnorm(M, mean = thetas[, 1], sd = 1)
Y <- rnorm(M, mean = thetas[, 2], sd = 1)


## PART 2 : Computation of the indices

# Build the matrix given by kernel evaluation and center them
K_XX <- do_outer_mat(X, X, "gaussian")
K_XX_center <- sweep(sweep(K_XX, 1, rowMeans(K_XX), "-"), 2, colMeans(K_XX), "-") + mean(K_XX)
K_YY <- do_outer_mat(Y, Y, "gaussian")
K_YY_center <- sweep(sweep(K_YY, 1, rowMeans(K_YY), "-"), 2, colMeans(K_YY), "-") + mean(K_YY)


## PART 2.1 : Pearson correlation coefficient
pear_corr <- cor(X, Y)


## PART 2.2 : Kernelized Canonical Correlation compute (KCC)

# regualirization parameter
epsilon <- 0.01

M1 <- chol2inv(chol(K_XX_center + epsilon * diag(n))) %*% K_XX_center
M2 <- chol2inv(chol(K_YY_center + epsilon * diag(n))) %*% K_YY_center

KCC <- svd(M1 %*% t(M2))$d[1]


# PART 2.3 : Constrained Covariance (COCO)

COCO <- sqrt(svd(K_XX_center %*% K_YY_center)$d[1]) / n


## PART 2.4 : Centered Kernel Alignment (CKA)

CKA <- sum(diag(K_XX_center %*% K_YY_center)) / sqrt(sum(diag(K_XX_center %*% K_XX_center)) * sum(diag(K_YY_center %*% K_YY_center)))


# store all the indices as a list
val_list <- list(pear_corr = pear_corr,
                 KCC = KCC,
                 COCO = COCO,
                 CKA = CKA)

# save the data matrix
write.csv(val_list, "output/results/other_indx_val_list.csv")
