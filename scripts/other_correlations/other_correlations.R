library(mvtnorm)

N <- 1000


rho_vec <- c(0,0.1,0.25,0.5,0.75,1)
rho_l <- length(rho_vec)


pear_corr_vec <- rep(NA,rho_l)

tr_corr_gau_vec <- rep(NA,rho_l)
CKA_gau_vec <- rep(NA,rho_l)

tr_corr_lin_vec <- rep(NA,rho_l)
CKA_lin_vec <- rep(NA,rho_l)


mat_center <- function(mat){
  mat <- sweep(sweep(mat, 1, rowMeans(mat), "-"), 2, colMeans(mat), "-") + mean(mat)
  return(mat)
}

tr_cov <- function(mat){
  n <- dim(mat)[1]
  return(tr(mat_center(mat))/(n-1))
}

HSIC <- function(mat1,mat2){
  n1 <- dim(mat1)[1]
  n2 <- dim(mat2)[1]
  return(tr(mat1 %*% mat2)/((n1-1)*(n2-1)))
}


for(i in seq_len(rho_l)){
  print(i)
  
  rho <- rho_vec[i]
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  
  set.seed(1234)
  thetas <- rmvnorm(N, sigma = sigma)
  
  X1 <- rnorm(N,mean = thetas[,1],sd = 1)
  X2 <- rnorm(N,mean = thetas[,1],sd = 1)
  X <- c(X1,X2)
  
  Y1 <- rnorm(N,mean = thetas[,2],sd = 1)
  Y2 <- rnorm(N,mean = thetas[,2],sd = 1)
  Y <- c(Y1,Y2)
  
  
  # Pearson correlation
  pear_corr_vec[i] <- cor(X,Y)
  
  
  # Gaussian kernel
  K_X12 <- do_outer_mat(X1,X2,"gaussian")
  K_Y12 <- do_outer_mat(Y1,Y2,"gaussian")
  
  K_XX <- do_outer_mat(X,X,"gaussian")
  K_YY <- do_outer_mat(Y,Y,"gaussian")
  K_XY <- do_outer_mat(X,Y,"gaussian")
  
  # Trace correlation
  tr_corr_gau_vec[i] <- tr_cov(K_XY)/sqrt(tr_cov(K_X12)*tr_cov(K_Y12))
  
  # CKA
  CKA_gau_vec[i] <- HSIC(K_XX,K_YY)/sqrt(HSIC(K_XX,K_XX)*HSIC(K_YY,K_YY))
  
  
  # Linear kernel
  K_X12 <- do_outer_mat(X1,X2,"linear")
  K_Y12 <- do_outer_mat(Y1,Y2,"linear")
  
  K_XX <- do_outer_mat(X,X,"linear")
  K_YY <- do_outer_mat(Y,Y,"linear")
  K_XY <- do_outer_mat(X,Y,"linear")
  
  # Trace correlation
  tr_corr_lin_vec[i] <- tr_cov(K_XY)/sqrt(tr_cov(K_X12)*tr_cov(K_Y12))
  
  # CKA
  CKA_gau_vec[i] <- HSIC(K_XX,K_YY)/sqrt(HSIC(K_XX,K_XX)*HSIC(K_YY,K_YY))
}


rho <- 1





# Pearson Correlation
pear_corr[i] <- cor(X,Y)

K_XX <- do_outer_mat(X,X,"gaussian")
K_XX_center <- sweep(sweep(K_XX, 1, rowMeans(K_XX), "-"), 2, colMeans(K_XX), "-") + mean(K_XX)
K_YY <- do_outer_mat(Y,Y,"gaussian")
K_YY_center <- sweep(sweep(K_YY, 1, rowMeans(K_YY), "-"), 2, colMeans(K_YY), "-") + mean(K_YY)


# HSIC
HSIC[i] <- tr(K_XX_center %*% K_YY_center)


# CKA
CKA[i] <- tr(K_XX_center %*% K_YY_center)/sqrt(tr(K_XX_center %*% K_XX_center)*tr(K_YY_center %*% K_YY_center))


# KCC
epsilon <- 0.01

M1 <- chol2inv(chol(K_XX_center + epsilon*diag(N))) %*% K_XX_center
M2 <- chol2inv(chol(K_YY_center + epsilon*diag(N))) %*% K_YY_center

KCC[i] <- svd(M1 %*% t(M2))$d[1]


# COCO

COCO[i] <- sqrt(svd(K_XX_center %*% K_YY_center)$d[1])/N


X <- rnorm(N, mean = thetas[,1], sd = 1)
Y <- rnorm(N, mean = thetas[,2], sd = 1)
rho_vec <- c(1)
rho_l <- length(rho_vec)

CKA <- rep(NA,10)
pear_corr <- rep(NA,10)
HSIC <- rep(NA,10)
CKA <- rep(NA,10)
KCC <- rep(NA,10)
COCO <- rep(NA,10)


for(i in seq_len(10)){
  print(i)
  rho <- 1
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  thetas <- rmvnorm(N, sigma = sigma)
  
  X <- rnorm(N, mean = thetas[,1], sd = 1)
  Y <- rnorm(N, mean = thetas[,2], sd = 1)

  
  # Pearson Correlation
  pear_corr[i] <- cor(X,Y)
  
  K_XX <- do_outer_mat(X,X,"gaussian")
  K_XX_center <- sweep(sweep(K_XX, 1, rowMeans(K_XX), "-"), 2, colMeans(K_XX), "-") + mean(K_XX)
  K_YY <- do_outer_mat(Y,Y,"gaussian")
  K_YY_center <- sweep(sweep(K_YY, 1, rowMeans(K_YY), "-"), 2, colMeans(K_YY), "-") + mean(K_YY)
  
  
  # HSIC
  HSIC[i] <- tr(K_XX_center %*% K_YY_center)
  
  
  # CKA
  CKA[i] <- tr(K_XX_center %*% K_YY_center)/sqrt(tr(K_XX_center %*% K_XX_center)*tr(K_YY_center %*% K_YY_center))
  
  
  # KCC
  epsilon <- 0.01
  
  M1 <- chol2inv(chol(K_XX_center + epsilon*diag(N))) %*% K_XX_center
  M2 <- chol2inv(chol(K_YY_center + epsilon*diag(N))) %*% K_YY_center
  
  KCC[i] <- svd(M1 %*% t(M2))$d[1]
  
  
  # COCO
  
  COCO[i] <- sqrt(svd(K_XX_center %*% K_YY_center)$d[1])/N
  

}

length(X)/2

X1 <- X[1:1000]
X2 <- X[1001:2000]
Y1 <- Y[1:1000]
Y2 <- Y[1001:2000]
K_XX <- do_outer_mat(X1,X2,"gaussian")
K_YY <- do_outer_mat(Y1,Y2,"gaussian")

K_XY <- do_outer_mat(X,Y,"gaussian")

cov <- tr(K_XY)/(2000 - 1) - sum(K_XY)/(2000*(2000 - 1))
varX <- tr(K_XX)/(1000 - 1) - sum(K_XX)/(1000*(1000 - 1))
varY <- tr(K_YY)/(1000 - 1) - sum(K_YY)/(1000*(1000 - 1))
cov/sqrt(varX*varY)

K_XX <- do_outer_mat(X1,X2,"gaussian")
K_XX_center <- sweep(sweep(K_XX, 1, rowMeans(K_XX), "-"), 2, colMeans(K_XX), "-") + mean(K_XX)
K_YY <- do_outer_mat(Y1,Y2,"gaussian")
K_YY_center <- sweep(sweep(K_YY, 1, rowMeans(K_YY), "-"), 2, colMeans(K_YY), "-") + mean(K_YY)
tr(K_XX_center %*% K_YY_center)/sqrt(tr(K_XX_center %*% K_XX_center)*tr(K_YY_center %*% K_YY_center))

mean(c(tr(K_XY1_center),tr(K_XY2_center)))/sqrt(tr(K_XX_center)*tr(K_YY_center))

N <- 1000

set.seed(1234)

rho_vec <- c(1)
rho_l <- length(rho_vec)

CA_k1 <- rep(NA,rho_l)
CA_k2 <- rep(NA,rho_l)
pear_corr1 <- rep(NA,rho_l)
pear_corr2 <- rep(NA,rho_l)

for(i in seq_len(rho_l)){
  print(i)
  rho <- rho_vec[i]
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  thetas <- rmvnorm(N, sigma = sigma)
  
  X1 <- rnorm(N, mean = thetas[,1], sd = 1)
  X2 <- rnorm(N, mean = thetas[,1], sd = 1)
  Y1 <- rnorm(N, mean = thetas[,2], sd = 1)
  Y2 <- rnorm(N, mean = thetas[,2], sd = 1)
  
  K_XX <- do_outer_mat(X1,X2,"gaussian")
  K_XX_center <- sweep(sweep(K_XX, 1, rowMeans(K_XX), "-"), 2, colMeans(K_XX), "-") + mean(K_XX)
  K_YY <- do_outer_mat(Y1,Y2,"gaussian")
  K_YY_center <- sweep(sweep(K_YY, 1, rowMeans(K_YY), "-"), 2, colMeans(K_YY), "-") + mean(K_YY)
  
  CA_k1[i] <- tr(K_XX_center %*% K_YY_center)/sqrt(tr(K_XX_center %*% K_XX_center)*tr(K_YY_center %*% K_YY_center))
  
  K_XX <- do_outer_mat(X1,X1,"gaussian")
  K_XX_center <- sweep(sweep(K_XX, 1, rowMeans(K_XX), "-"), 2, colMeans(K_XX), "-") + mean(K_XX)
  K_YY <- do_outer_mat(Y1,Y1,"gaussian")
  K_YY_center <- sweep(sweep(K_YY, 1, rowMeans(K_YY), "-"), 2, colMeans(K_YY), "-") + mean(K_YY)
  
  CA_k2[i] <- tr(K_XX_center %*% K_YY_center)/sqrt(tr(K_XX_center %*% K_XX_center)*tr(K_YY_center %*% K_YY_center))
  
  pear_corr1[i] <- cor(X1,Y1)
  pear_corr2[i] <- cov(X1,Y1)/sqrt(cov(X1,X2)*cov(Y1,Y2))
  
}