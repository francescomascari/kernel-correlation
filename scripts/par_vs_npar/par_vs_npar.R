## Sampling functions:

hDP_mat_X2_help <- function(
  n,
  group_indx,
  seen = matrix(numeric(0),ncol=4,dimnames = list(NULL,c("X", "Ncusts1", "Ncusts2", "Ntabs"))),
  c0 = 1,
  c = 1,
  P00 = runif){
  # ---------------------------------------------------------------------------
  # Draw n observations for a single group from the hierarchical DP (HDP) model.
  # Sampling follows the Restaurant franchise scheme with table (cluster) indicators:
  #   * A new customer joins an existing table with prob. Ncusts / (c + Ncusts)
  #   * Otherwise a new table is created.
  #       – The table serves an existing dish with prob. Ntabs / (c0 + Ntabs)
  #       – Otherwise it introduces a brand-new dish from the baseline measure P00.
  # The function updates the shared `seen` count matrix and returns the updated
  # matrix plus the n newly generated (value, group) rows.
  # ---------------------------------------------------------------------------
  #
  # Arguments
  #   n           : integer – number of customers to add to the given group.
  #
  #   group_indx  : integer (1 or 2) – which group these customers belong to.
  #
  #   seen        : numeric matrix with columns
  #                   X        – unique dish values
  #                   Ncusts1  – frequency in group 1
  #                   Ncusts2  – frequency in group 2
  #                   Ntabs    – number of tables serving the dish
  #                 If empty, the process starts from scratch.
  #
  #   c, c0       : numeric – concentration parameters of the HDP.
  #
  #   P00         : function – baseline measure generator (default `runif`).
  #
  # Returns
  #   list(
  #     seen   = updated `seen` matrix (same column structure as above),
  #     new_X  = matrix of newly sampled rows with columns:S
  #                X     – sampled value
  #                group – group index (matches `group_indx`)
  #   )
  
  # Initialize the matrix of newly sampled values
  new_X <-  matrix(numeric(0),ncol=2,dimnames = list(NULL,c("X", "group")))
  
  # If no value is requested to be sampled,
  # return the original `seen` matrix and the empty `new_X` matrix
  if(n==0){
    return(list(seen=seen,new_X=new_X))
  }
  
  # For the data in the `seen` matrix, record
  n_dishes <- nrow(seen) #number of unique values already sampled
  n_tabs <- sum(seen[,"Ntabs"]) #number of tables
  n_custs <- sum(seen[,group_indx+1]) #number of costumers
  
  # For every new value to be sampled
  for (k in seq_len(n)) {
    if (runif(1) < c/(c + n_custs)) {
      # new table
      if (runif(1) < c0/(c0 + n_tabs)){
        # new dish
        dish <- P00(1)
        n_dishes <- n_dishes + 1
        
        # add the entry related to the new dish to the `seen` matrix
        new_dish <- c(dish,0,0,1)
        new_dish[group_indx+1] <- 1
        seen <- rbind(seen,new_dish)
      }
      else{
        # old dish
        
        # sample an index in seen with respect to the distribution of tables
        indx <- sample(1:n_dishes,size = 1,prob = seen[,"Ntabs"])
        
        # update the entry related to the sampled dish in the `seen` matrix
        dish <- seen[indx,"X"]
        seen[indx,"Ntabs"] <- seen[indx,"Ntabs"] + 1
        seen[indx,group_indx+1] <- seen[indx,group_indx+1] + 1
      }
      n_tabs <- n_tabs + 1
    }
    else {
      # old table, old dish
      
      # sample an index in seen with respect to the distribution of customers in the group
      indx <- sample(1:n_dishes,size = 1,prob = seen[,group_indx+1])
      
      # update the entry related to the sampled dish in the `seen` matrix
      dish <- seen[indx,"X"]
      seen[indx,group_indx+1] <- seen[indx,group_indx+1] + 1
    }
    n_custs <- n_custs + 1
    
    # update the `new_X` matrix with the sampled value and the corresponding group
    new_X <- rbind(new_X,c(dish,group_indx))
  }
  
  return(list(seen = seen,new_X = new_X))
}


hDP_X2_mat_data <- function(
  n1,
  n2,
  smpl_method = "block",
  start = 1,
  seen = matrix(numeric(0),ncol=4,dimnames = list(NULL,c("X", "Ncusts1", "Ncusts2", "Ntabs"))),
  c = 1,
  c0 = 1,
  P00 = runif){
  # ---------------------------------------------------------------------------
  # Sample n1 and n2 observations from a two-group hierarchical DP (HDP) model.
  # If `seen` is supplied, draws are conditional on the counts in that matrix.
  # The augmented model with tables is used.
  # ---------------------------------------------------------------------------
  #
  # Arguments
  #   n1, n2      : integers – numbers of observations to draw from groups 1 and 2.
  #
  #   smpl_method : character – sampling schedule
  #                 * "block" – sample all from one group, then all from the other,
  #                             starting with the group indicated by `start`.
  #                 * "alt"   – alternate single observations, starting with `start`.
  #
  #   start       : integer (1 or 2) – group to draw the first observation from.
  #
  #   seen        : numeric matrix with columns
  #                   X        – unique sampled values
  #                   Ncusts1  – frequency in group 1
  #                   Ncusts2  – frequency in group 2
  #                   Ntabs    – number of CRP tables
  #                 If empty, sampling is unconditional.
  #
  #   c, c0       : numeric – concentration parameters of the HDP.
  #
  #   P00         : function – baseline measure generator (default `runif`).
  #
  # Returns
  #   list(
  #     seen   = updated `seen` matrix (same column structure as above),
  #     new_X  = matrix of newly sampled rows with columns:
  #                X     – sampled value
  #                group – group index (1 or 2)
  #   )
  
  # Assign the first and the second indices and the corresponding cardinality
  # Based on the value of `start`.
  if (start == 1) {
    indx_1st <- 1
    n_1st <- n1
    
    indx_2nd <- 2
    n_2nd <- n2
  }
  else {
    indx_1st <- 2
    n_1st <- n2
    
    indx_2nd <- 1
    n_2nd <- n1
  }
  
  # Initialize the matrix of newly sampled values
  new_X <- matrix(numeric(0),ncol=2,dimnames = list(NULL,c("X", "group")))
  
  if (smpl_method == "block") {
    # If the sampling method is "block"
    # we sample all values from the first group
    group1_updt <- hDP_mat_X2_help(n_1st,group_indx = indx_1st, seen = seen, c0 = c0, c = c,P00 = P00)
    seen <- group1_updt$seen
    new_X <- rbind(new_X,group1_updt$new_X)
    
    # we sample all values from the second group
    group2_updt <- hDP_mat_X2_help(n_2nd,group_indx = indx_2nd, seen = seen, c0 = c0, c = c,P00 = P00)
    seen <- group2_updt$seen
    new_X <- rbind(new_X,group2_updt$new_X)
  }
  else if (smpl_method == "alt") {
    # If the sampling method is "block"
    
    n_min <- min(n1,n2) # minimal number of values to sample
    n_max <- max(n1,n2) # maximal number of values to sample
    max_indx <- which.max(c(n1,n2)) # group corresponding to the maximal number of values to samples
    
    #until we have to sample from both groups
    for (k in seq_len(n_min)) {
      # we sample one value from the first group
      group1_updt <- hDP_mat_X2_help(1,group_indx = indx_1st, seen = seen, c0 = c0, c = c,P00 = P00)
      seen <- group1_updt$seen
      new_X <- rbind(new_X,group1_updt$new_X)
      
      # we sample one value from the second group
      group2_updt <- hDP_mat_X2_help(1,group_indx = indx_2nd, seen = seen, c0 = c0, c = c,P00 = P00)
      seen <- group2_updt$seen
      new_X <- rbind(new_X,group2_updt$new_X)
    }
    
    # we sample the remaining ones from the group corresponding to the maximal number of values
    last_updt <- hDP_mat_X2_help(n_max - n_min,group_indx = max_indx, seen = seen, c0 = c0, c = c,P00 = P00)
    seen <- last_updt$seen
    new_X <- rbind(new_X,last_updt$new_X)
  }
  return(list(seen=seen,new_X=new_X))
}


## Posterior analysis

compute_s_sq <- function(t_sq, v, sigma){
  return(sigma^2*((v + sqrt(sigma^2/(2*t_sq + sigma^2)))^-2 - 1)/2)
}

# compute_gamma_sq <- function(t_sq, s_sq){
#   return(t_sq - s_sq)
# }

compute_gamma_sq <- function(t_sq, v, sigma){
  return(t_sq - sigma^2*((v + sqrt(sigma^2/(2*t_sq + sigma^2)))^-2 - 1)/2)
}

# compute_rho <- function(xi, gamma_sq, s_sq, sigma){
#   # inner term
#   inner <- xi + (1 - xi)*sqrt((2*s_sq + sigma^2)/(2*gamma_sq + 2*s_sq + sigma^2))
# 
#   # apply the formula
#   rho <- 1 - (2*s_sq + sigma^2)/(2*gamma_sq)*(inner^(-2) - 1)
# 
#   return(rho)
# }

compute_rho <- function(xi, t_sq, v, sigma){
  # inner term
  up <- t_sq - sigma^2*((v*xi + sqrt(sigma^2/(2*t_sq + sigma^2)))^-2 - 1)/2
  
  # apply the formula
  down <- t_sq - sigma^2*((v + sqrt(sigma^2/(2*t_sq + sigma^2)))^-2 - 1)/2
  
  return(up/down)
}

# compute_c0 <- function(gamma_sq, rho, s_sq, sigma){
#   # intermediate terms
#   t1 <- sqrt(sigma^2/(2*gamma_sq*(1 - rho) + 2*s_sq + sigma^2))
#   t2 <- sqrt(sigma^2/(2*gamma_sq + 2*s_sq + sigma^2))
#   # assemble c0
#   return((1 - t1)/(t1 - t2))
# }

compute_c0 <- function(xi, t_sq, v, sigma){
  return((1 - sqrt(sigma^2/(2*t_sq + sigma^2)))/(v*xi) - 1)
}

# compute_c <- function(gamma_sq, rho, s_sq, sigma){
#   # intermediate terms
#   t1 <- sqrt(sigma^2/(2*s_sq + sigma^2))
#   t2 <- sqrt(sigma^2/(2*gamma_sq*(1 - rho) + 2*s_sq + sigma^2))
#   
#   # numerator and denominator
#   return((1 - t1)/(t1 - t2))
# }

compute_c <- function(xi, t_sq, v, sigma){
  return(((1 - sqrt(sigma^2/(2*t_sq + sigma^2)))/v - 1)/(1-xi))
}

## PLOT FOR RHO, c0 AND C

t_sq <- 2
v <- 1/4

sigma_max <- sqrt(2*t_sq/(1/(1-v)^2 - 1))
sigma_vec <- sigma_max*sqrt(2)^seq(-1,-20,by = -4)
xi_vec <- c(0.01,seq(0,1,length.out = 21)[1:21],0.99)

### PLOT FOR RHO


rho_mat <- sapply(sigma_vec, function(sigma) compute_rho(xi_vec, t_sq, v, sigma), USE.NAMES = TRUE)
rho_df <- as.data.frame(rho_mat)
names(rho_df) <- c(as.character(format(round(sigma_vec,2),nsmall = 2)))
rho_df$corr <- xi_vec

rho_df_long <- rho_df %>%
  pivot_longer(cols = -corr, 
               names_to = "sigma", 
               values_to = "rho") %>%
  mutate(sigma = as.factor(sigma))

ggplot(rho_df_long, aes(x = corr, y = rho, color = sigma, shape = sigma)) +
  geom_line(linewidth = 2) +
  geom_point(size = 6) +
  labs(y = expression(rho), 
       x = "Correlation") +
  #scale_color_manual(values = c("0.01"="purple3", "0.1"="darkorange3","1"="lightgreen","10"="cyan3","100"="darkred")) + # associate each value of fill to a color
  #scale_shape_manual(values = c("0.01"=15, "0.1"=16,"1"=18,"10"=17,"100"=20)) +
  theme_classic() + 
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
        panel.grid = element_line(size = 1.5),
        legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(family="LM Roman 10",size=20,face="bold"))

## PLOT FOR c0

c0_mat <- sapply(sigma_vec, function(sigma) compute_c0(xi_vec, t_sq, v, sigma), USE.NAMES = TRUE)
c0_df <- as.data.frame(c0_mat)
names(c0_df) <- c(as.character(format(round(sigma_vec,2),nsmall = 2)))
c0_df$corr <- xi_vec

c0_df_long <- c0_df %>%
  pivot_longer(cols = -corr, 
               names_to = "sigma", 
               values_to = "c0") %>%
  mutate(sigma = as.factor(sigma))

ggplot(c0_df_long, aes(x = corr, y = c0, color = sigma, shape = sigma)) +
  geom_line(linewidth = 2) +
  geom_point(size = 6) +
  labs(y = expression(c[0]), 
       x = "Correlation") +
  scale_y_continuous(transform = "log10") +
  #scale_color_manual(values = c("0.01"="purple3", "0.1"="darkorange3","1"="lightgreen","10"="cyan3","100"="darkred")) + # associate each value of fill to a color
  #scale_shape_manual(values = c("0.01"=15, "0.1"=16,"1"=18,"10"=17,"100"=20)) +
  theme_classic() + 
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
        panel.grid = element_line(size = 1.5),
        legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(family="LM Roman 10",size=20,face="bold"))

## PLOT FOR c

c_mat <- sapply(sigma_vec, function(sigma) compute_c(xi_vec, t_sq, v, sigma), USE.NAMES = TRUE)
c_df <- as.data.frame(c_mat)
names(c_df) <- c(as.character(format(round(sigma_vec,2),nsmall = 2)))
c_df$corr <- xi_vec

c_df_long <- c_df %>%
  pivot_longer(cols = -corr, 
               names_to = "sigma", 
               values_to = "c") %>%
  mutate(sigma = as.factor(sigma))

ggplot(c_df_long, aes(x = corr, y = c, color = sigma, shape = sigma)) +
  geom_line(linewidth = 2) +
  geom_point(size = 6) +
  labs(y = expression(c), 
       x = "Correlation") +
  scale_y_continuous(transform = "log10") +
  #scale_color_manual(values = c("0.01"="purple3", "0.1"="darkorange3","1"="lightgreen","10"="cyan3","100"="darkred")) + # associate each value of fill to a color
  #scale_shape_manual(values = c("0.01"=15, "0.1"=16,"1"=18,"10"=17,"100"=20)) +
  theme_classic() + 
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
        panel.grid = element_line(size = 1.5),
        legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(family="LM Roman 10",size=20,face="bold"))

### DATA GENERATION

t_sq <- 2
v <- 1/4
sigma <- sqrt(t_sq/(1/(1-v)^2 - 1))

s_sq <- compute_s_sq(t_sq, v, sigma)
gamma_sq <- compute_gamma_sq(t_sq, v, sigma)


set.seed(1)
n1 <- 200
X1 <- hDP_XT2(n1,0,c=10,c0=10,P00=function(n){rnorm(n,mean = -1, sd = sqrt(s^{2} + g^{2}))})[,"X"]
m1 <- mean(X1)

n2 <- 10
X2 <- hDP_XT2(0,n2,c=10,c0=10,P00=function(n){rnorm(n,mean = 1, sd = sqrt(s^{2} + g^{2}))})[,"X"]
m2 <- mean(X2)

X_tot <- unique(c(X1,X2))
Ncusts1 <- sapply(X_tot,function(val){sum(X1 == val)})
Ncusts2 <- sapply(X_tot,function(val){sum(X2 == val)})
seen <- data.frame(X=X_tot,Ncusts1=Ncusts1,Ncusts2=Ncusts2)

M <- 10000

xi_vec <- c(0.001,0.50,0.99)
xi_len <- length(xi_vec)

## Gaussian

X_gau <- data.frame(X = numeric(0),
                    Group = factor(levels = c("Group 1", "Group 2")), 
                    Corr = factor(levels = as.character(xi_vec))
)

for(i in seq_len(xi_len)){
  corr <- xi_vec[i]
  rho <- compute_rho(corr, t_sq, v, sigma)

  #Compute the marginal means
  mu1 <- (n1*m1 + n1*n2*gamma_sq/s_sq*(1-rho^2)*m1 + rho*n2*m2)/(s_sq/gamma_sq + n1 + n2 + n1*n2*gamma_sq/s_sq*(1-rho^2))
  mu2 <- (n2*m2 + n1*n2*gamma_sq/s_sq*(1-rho^2)*m2 + rho*n1*m1)/(s_sq/gamma_sq + n1 + n2 + n1*n2*gamma_sq/s_sq*(1-rho^2))
  
  #Compute the marginal variances
  V1 <- s_sq*(s_sq/gamma_sq + n1 + n2 + 1 + (n1+1)*n2*gamma_sq/s_sq*(1-rho^2))/(s_sq/gamma_sq + n1 + n2 + n1*n2*gamma_sq/s_sq*(1-rho^2))
  V2 <- s_sq*(s_sq/gamma_sq + n1 + n2 + 1 + n1*(n2+1)*gamma_sq/s_sq*(1-rho^2))/(s_sq/gamma_sq + n1 + n2 + n1*n2*gamma_sq/s_sq*(1-rho^2))
  
  set.seed(2)
  X1_vals <- rnorm(M,mean = mu1, sd = sqrt(V1))
  X2_vals <- rnorm(M,mean = mu2, sd = sqrt(V2))
  
  X_gau <- rbind(X_gau,
                 data.frame(X = X1_vals,
                            Group = factor("Group 1",levels = c("Group 1", "Group 2")),
                            Corr = factor(as.character(corr),levels = as.character(xi_vec))),
                 data.frame(X = X2_vals,
                            Group = factor("Group 2",levels = c("Group 1", "Group 2")),
                            Corr = factor(as.character(corr),levels = as.character(xi_vec)))
  )
}


## hDP

X_vec <- seen[,1]

n1_vec <- seen[,2]
n1_tot <- sum(n1_vec)

n2_vec <- seen[,3]
n2_tot <- sum(n2_vec)

if(n1_tot != 0 || n2_tot != 0){
  seen_q1_mat <- lapply(n1_vec, FUN = function(times){rep(1,times)})
  seen_q2_mat <- lapply(n2_vec, FUN = function(times){rep(1,times)})
}

X_hDP <- data.frame(X = numeric(0),
                    Group = factor(levels = c("Group 1", "Group 2")), 
                    Corr = factor(levels = as.character(xi_vec))
)

for(i in seq_len(xi_len)){
  corr <- xi_vec[i]
  c0 <- compute_c0(corr, t_sq, v, sigma)
  c <- compute_c(corr, t_sq, v, sigma)
  
  set.seed(2)
  
  X1_vals <- c()
  X2_vals <- c()
  
  for(i in seq_len(M)){
    l_vec <- sapply(seen_q1_mat,length) + sapply(seen_q2_mat,length)
    seen_mat <- cbind(seen,"Ntabs"=l_vec) # we impose this table disposition
    
    joint_vals1 <- hDP_X2_mat_data(1,0,seen = seen_mat,smpl_method = "alt",start = 1,c=c,c0=c0,P00=function(n){rnorm(n,mean = 0, sd = sqrt(t_sq))})$new_X
    X1_vals <- c(X1_vals,joint_vals1[joint_vals1[,"group"]==1,1][1])
    joint_vals2 <- hDP_X2_mat_data(0,1,seen = seen_mat,smpl_method = "alt",start = 2,c=c,c0=c0,P00=function(n){rnorm(n,mean = 0, sd = sqrt(t_sq))})$new_X
    X2_vals <- c(X2_vals,joint_vals2[joint_vals2[,"group"]==2,1][1])
    
    seen_q1_mat <- gibbs_tabs(l_vec,seen_q1_mat,c0=c0,c=c)
    l_vec <- sapply(seen_q1_mat,length) + sapply(seen_q2_mat,length)
    seen_q2_mat <- gibbs_tabs(l_vec,seen_q2_mat,c0=c0,c=c)
  }
  
  X_hDP <- rbind(X_hDP,
                 data.frame(X = X1_vals,
                            Group = factor("Group 1",levels = c("Group 1", "Group 2")),
                            Corr = factor(as.character(corr),levels = as.character(xi_vec)))
  )
  X_hDP <- rbind(X_hDP,
                 data.frame(X = X2_vals,
                            Group = factor("Group 2",levels = c("Group 1", "Group 2")),
                            Corr = factor(as.character(corr),levels = as.character(xi_vec)))
  )
  
}

#Plots and

for(i in seq_len(xi_len)){
  corr <- xi_vec[i]
  
  a <- ggplot(subset(X_gau,Corr == as.character(corr)), aes(x = X,color = Group, linetype = Group)) +
    geom_density(position = "identity", fill = scales::alpha("darkgray", 0.66), linewidth = 2) +
    geom_density(position = "identity", fill = NA, linewidth = 2) +
    scale_color_manual(values = c("Group 1"="firebrick", "Group 2"="forestgreen")) +
    scale_linetype_manual(values = c("Group 1"="solid","Group 2"="dashed")) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
          panel.grid = element_line(size = 1.5),
          legend.position="top",
          legend.title = element_blank(),
          legend.text = element_text(family="LM Roman 10",size=20,face="bold")) +
    xlim(c(-4,4)) +
    ylim(c(0,2))
  
  b <- ggplot(subset(X_hDP,Corr == as.character(corr)), aes(x = X,y = after_stat(density),color = Group, linetype = Group)) +
    geom_histogram(position = "identity", fill = scales::alpha("darkgray", 0.66), linewidth = 2) + 
    geom_histogram(position = "identity", fill = NA, linewidth = 2) +
    scale_color_manual(values = c("Group 1"="firebrick", "Group 2"="forestgreen")) +
    scale_linetype_manual(values = c("Group 1"="solid","Group 2"="dashed")) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
          panel.grid = element_line(size = 1.5),
          legend.position="top",
          legend.title = element_blank(),
          legend.text = element_text(family="LM Roman 10",size=20,face="bold")) +
    ylim(c(0,1.2)) +
    xlim(c(-5.5,5.5))
  
  show(a)
  show(b)
}

corr_df <- expand.grid(corr1 = as.factor(xi_vec), corr2 = as.factor(xi_vec))

val_gau <- numeric(xi_len^2)
val_hDP <- numeric(xi_len^2)
val_gauVShDP <- numeric(xi_len^2)

for(i in seq_len(xi_len^2)){
  corr1 <- corr_df$corr1[i]
  corr2 <- corr_df$corr2[i]
  
  val_gau[i] <- abs(mean(subset(X_gau,Corr == as.character(corr1) & Group == "Group 2")$X) - mean(subset(X_gau,Corr == as.character(corr2) & Group == "Group 2")$X))
  val_hDP[i] <- abs(mean(subset(X_hDP,Corr == as.character(corr1) & Group == "Group 2")$X) - mean(subset(X_hDP,Corr == as.character(corr2) & Group == "Group 2")$X))
  val_gauVShDP[i] <- abs(mean(subset(X_gau,Corr == as.character(corr1) & Group == "Group 2")$X) - mean(subset(X_hDP,Corr == as.character(corr2) & Group == "Group 2")$X))
}

gau_df <- cbind(corr_df,value = val_gau)
hDP_df <- cbind(corr_df,value = val_hDP)
gauVShDP_df <- cbind(corr_df,value = val_gauVShDP)


a <- ggplot(gau_df, aes(x = corr1, y = corr2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", family="LM Roman 10", size = 8, fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient(high="firebrick", low="forestgreen") +
  labs(x = "Correlation: Gaussian case", 
       y = "Correlation: Gaussian case") +
  theme_classic() + 
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
        panel.grid = element_line(size = 1.5),
        legend.position="none")

b <- ggplot(hDP_df, aes(x = corr1, y = corr2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", family="LM Roman 10", size = 8, fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient(high="firebrick", low="forestgreen") +
  labs(x = "Correlation: hDP case", 
       y = "Correlation: hDP case") +
  theme_classic() + 
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
        panel.grid = element_line(size = 1.5),
        legend.position="none")

c <- ggplot(gauVShDP_df, aes(x = corr1, y = corr2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", family="LM Roman 10", size = 8, fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient(high="firebrick", low="forestgreen") +
  labs(x = "Correlation: Gaussian case", 
       y = "Correlation: hDP case") +
  theme_classic() + 
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
        panel.grid = element_line(size = 1.5),
        legend.position="none")
  
show(a)
show(b)
show(c)



###### EXTRA ######

## Prior analysis hDP

sigma <-1 
t_sq <- 2
v <- 1/2
s_sq <- compute_s_sq(t_sq, v, sigma)
gamma_sq <- compute_gamma_sq(t_sq, s_sq)

corr_vec <- c(0.5)
for(corr in corr_vec){
  rho <- compute_rho(corr, gamma_sq, s_sq, sigma)
  c0 <- compute_c0(gamma_sq, rho, s_sq, sigma)
  c <- compute_c(gamma_sq, rho, s_sq, sigma)
  
  M <- 10000
  X11_vals <- c()
  X12_vals <- c()
  X21_vals <- c()
  X22_vals <- c()
  for(i in seq_len(M)){
    joint_vals <- hDP_X2_mat_data(2,2,smpl_method = "alt",start = 1,c=c,c0=c0,P00=function(n){rnorm(n,mean = 0, sd = sqrt(gamma_sq + s_sq))})$new_X
    X11_vals <- c(X11_vals,joint_vals[joint_vals[,"group"]==1,1][1])
    X12_vals <- c(X12_vals,joint_vals[joint_vals[,"group"]==1,1][2])
    X21_vals <- c(X21_vals,joint_vals[joint_vals[,"group"]==2,1][1])
    X22_vals <- c(X22_vals,joint_vals[joint_vals[,"group"]==2,1][2])
  }
  
  true_varA <- sqrt(sigma^2 / (2*s_sq + sigma^2)) *(1 - sqrt((2*s_sq + sigma^2) / (2*gamma_sq + 2*s_sq + sigma^2)))
  true_varB <- ((1 + c + c0) / ((1 + c) * (1 + c0))) *(1 - sqrt(sigma^2 / (2*gamma_sq + 2*s_sq + sigma^2)))
  smpl_var1 <- (sum(do_outer_diag(X11_vals,X12_vals)) - sum(do_outer_mat(X11_vals,X12_vals))/M)/(M-1)
  smpl_var2 <- (sum(do_outer_diag(X21_vals,X22_vals)) - sum(do_outer_mat(X21_vals,X22_vals))/M)/(M-1)
  
  print(c(true_varA,true_varB,smpl_var1))
  print(c(true_varA,true_varB,smpl_var2))
  
  true_covA <- sqrt(sigma^2 / (2*s_sq + sigma^2)) *(sqrt((2*s_sq + sigma^2) / (2*gamma_sq * (1 - rho) + 2*s_sq + sigma^2)) - sqrt((2*s_sq + sigma^2) / (2*gamma_sq + 2*s_sq + sigma^2)))
  true_covB <- (1 / (1 + c0)) * (1 - sqrt(sigma^2 / (2*gamma_sq + 2*s_sq + sigma^2)))
  smpl_cov <- (sum(do_outer_diag(X11_vals,X21_vals)) - sum(do_outer_mat(X11_vals,X21_vals))/M)/(M-1)
  
  print(c(true_covA,true_covB,smpl_cov))
  
  print(c(corr,true_covA/true_varA,true_covB/true_varB,smpl_cov/sqrt(smpl_var1*smpl_var2))) 
}

## GAUSSIAN PRIOR

corr_gau <- function(rho,sigma = 1){
  return ((sqrt((4+sigma^2)/(4+sigma^2-2*rho))-1)/(sqrt((4+sigma^2)/(2+sigma^2))-1))
}

sigma_vec <- as.character(10^seq(-2,2,by = 1))
rho_vec <- seq(-1,1,length.out = 21)

corr_gau_mat <- sapply(sigma_vec, function(s) corr_gau(rho_vec, as.numeric(s)), USE.NAMES = TRUE)
corr_gau_df <- as.data.frame(corr_gau_mat)
corr_gau_df$rho <- rho_vec

corr_gau_df_long <- corr_gau_df %>%
  pivot_longer(cols = -rho, 
               names_to = "sigma", 
               values_to = "value")

corr_gau_df_long <- corr_gau_df_long %>%
  mutate(sigma = as.numeric(sigma))

ggplot(corr_gau_df_long, aes(x = rho, y = value, color = as.factor(sigma), shape = as.factor(sigma))) +
  geom_line(linewidth = 2) +
  geom_point(size = 6) +
  labs(x = expression(rho), 
       y = "Correlation") +
  scale_color_manual(values = c("0.01"="purple3", "0.1"="darkorange3","1"="lightgreen","10"="cyan3","100"="darkred")) + # associate each value of fill to a color
  scale_shape_manual(values = c("0.01"=15, "0.1"=16,"1"=18,"10"=17,"100"=20)) +
  theme_classic() + 
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
        panel.grid = element_line(size = 1.5),
        legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(family="LM Roman 10",size=20,face="bold"))



## FARE GEAFICI CON CORRELAZIONE SULLE X
## E VALORI DEI PARAMETRI PER OTTENERE TALE CORR
## SULLE Y
## GAU: CORR VS RHO
## HDP: CORR VS C0 (C FISSATO)
## HDP: CORR VS C (C0 FISSATO)

## hDP PRIOR

corr_hDP <- function(c0,c){
  return ((1 + c)/(1 + c + c0))
}

c0c_vec <- 10^seq(-2,2,by = 1)

corr_hDP_val <- outer(c0c_vec, c0c_vec, FUN = function(c0, c) corr_hDP(c0, c))

corr_hDP_df <- expand.grid(c0 = c0c_vec, c = c0c_vec)
corr_hDP_df$value <- as.vector(corr_hDP_val)


ggplot(corr_hDP_df, aes(x = as.factor(c0), y = as.factor(c), fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", family="LM Roman 10", size = 6, fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient(low="firebrick", high="forestgreen") +
  labs(x = bquote(bold(c[0])), 
       y = "c") +
  theme_classic() + 
  theme(axis.title = element_text(vjust = 0, family="LM Roman 10", size = 30,face="bold"),
        axis.text = element_text(size = 20, family="LM Roman 10", face = "bold"),
        panel.grid = element_line(size = 1.5),
        legend.position="none")

