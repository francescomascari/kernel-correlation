hDP_corr_analytics <- function(
  seen = matrix(numeric(0),ncol=3,dimnames = list(NULL,c("X", "Ncusts1", "Ncusts2"))),
  c0 = 1,
  c = 1,
  R = 1000,
  baseln = runif,
  M = 10000,
  ...){
  # ---------------------------------------------------------------------------
  # Compute the kernel correlation for the hierarchical DP (HDP) model using
  # using the "analytics" method specified in the supplementary material
  # of the manuscript.
  # If `seen` is supplied, the computation is a posteriori,
  # conditionally on the counts in that matrix.
  # ---------------------------------------------------------------------------
  #
  # Arguments
  #   seen        : numeric matrix with columns
  #                 (using the restaurant franchise metaphor)
  #                   X        – unique dish values
  #                   Ncusts1  – frequency of customers in group 1
  #                   Ncusts2  – frequency of customers in group 2
  #                 If empty (default), the computation is a priori.
  #
  #   c0, c       : numeric – concentration parameters of the HDP (default 1).
  #
  #   R           : integer – number of runs of the Gibbs sampler to
  #                           integrate out the tables.
  #
  #   baseln      : function – baseline measure generator (default `runif`).
  #
  #   M           : integer – size of the sample from the baseline measure.
  #
  #   ...         : extra arguments: kernel and its parameters
  #
  # Returns:
  #   numeric – the value of the kernel correlation.

    # sample M values from `baseln`
    baseln_smpl <- baseln(M)
    
    # store the unique values of the dishes
    X_vec <- seen[,"X"]
    
    # store the frequencies of the unique values for the first group
    n1_vec <- seen[,"Ncusts1"]
    n1_tot <- sum(n1_vec) # total number of customers for the first group
    
    # store the frequencies of the unique values for the second group
    n2_vec <- seen[,"Ncusts2"]
    n2_tot <- sum(n2_vec) # total number of customers for the second group
    
    if (n1_tot != 0 || n2_tot != 0) { # if at least one group has some values

      # assign one table per each customer for each group
      seen_q1_mat <- lapply(n1_vec, FUN = function(times){rep(1,times)})
      seen_q2_mat <- lapply(n2_vec, FUN = function(times){rep(1,times)})
      
      # pre-compute the kernel matrices for the computation of the integrals
      outer_mat_baseln <- do_outer_mat(vals1 = baseln_smpl,vals2 = baseln_smpl,...)
      outer_mat_dishes <- do_outer_mat(vals1 = X_vec,vals2 = X_vec,...)
      outer_mat_baseln_dishes <- do_outer_mat(vals1 = baseln_smpl,vals2 = X_vec,...)
      
      # pre-compute the kernel vectors for the computation of the integrals
      diag_vec_baseln <- do_diag_vec(vals = baseln_smpl,...)
      diag_vec_dishes <- do_diag_vec(vals = X_vec,...)
      
      # compute the double integral of the kernel k(x,y) w.r.t.
      # the baseline measure both in x and in y
      int_cross_baseln <- sum(outer_mat_baseln)/M^2
      
      # compute the integral of the kernel k(x,x) w.r.t
      # the baseline measure in x
      int_diag_baseln <- sum(diag_vec_baseln)/M
      
      if (n1_tot != 0) { # if at least one customer in group 1

        # compute the integral of the kernel k(x,y) w.r.t.
        # the empirical distribution of the dishes
        # according to the frequency of the customers in group 1 both in x and in y
        int_cross_custs1 <- quad.form(outer_mat_dishes,n1_vec)/n1_tot^2

        # compute the integral of the kernel k(x,y) w.r.t.
        # the baseline measure in x,
        # the empirical distribution of the dishes
        # according to the frequency of the customers in group 1 in y
        int_cross_baseln_custs1 <- sum(colSums(outer_mat_baseln_dishes)*n1_vec)/(M*n1_tot)
        
        # compute the integral of the kernel k(x,x) w.r.t
        # the empirical distribution of the dishes
        # according to the frequency of the customers in group 1 in x
        int_diag_custs1 <- sum(diag_vec_dishes*n1_vec)/n1_tot
      }
      else{ # if no customers in group 1

        # all the integrals involving the empirical distribution of the dishes
        # according to the frequency of the customers in group 1 are null
        int_cross_custs1 <- 0
        int_cross_baseln_custs1 <- 0
        int_diag_custs1 <- 0
      }
      
      # third component of the variance for group 1
      V1_3 <- n1^2/((c + n1 + 1)*(c + n1)^2)*(int_diag_custs1 - int_cross_custs1)
      
      if (n2_tot != 0) { # if at least one customer in group 2

        # compute the integral of the kernel k(x,y) w.r.t.
        # the empirical distribution of the dishes
        # according to the frequency of the customers in group 2 both in x and in y
        int_cross_custs2 <- quad.form(outer_mat_dishes,n2_vec)/n2_tot^2

        # compute the integral of the kernel k(x,y) w.r.t.
        # the baseline measure in x,
        # the empirical distribution of the dishes
        # according to the frequency of the customers in group 2 in y
        int_cross_baseln_custs2 <- sum(colSums(outer_mat_baseln_dishes)*n2_vec)/(M*n2_tot)
        
        # compute the integral of the kernel k(x,x) w.r.t
        # the empirical distribution of the dishes
        # according to the frequency of the customers in group 2 in x
        int_diag_custs2 <- sum(diag_vec_dishes*n2_vec)/n2_tot
      }
      else{ # if no customers in group 2

        # all the integrals involving the empirical distribution of the dishes
        # according to the frequency of the customers in group 2 are null
        int_cross_custs2 <- 0
        int_cross_baseln_custs2 <- 0
        int_diag_custs2 <- 0
      }

      # third component of the variance for group 2
      V2_3 <- n2^2/((c + n2 + 1)*(c + n2)^2)*(int_diag_custs2 - int_cross_custs2)
      
      # pre-allocate the values for which we have to integrate out the tables
      V0_1 <- 0
      V1_1 <- 0
      V1_2 <- 0
      V2_1 <- 0
      V2_2 <- 0
      
      for (. in seq_len(R)) { # Gibbs iterations

        # store the frequencies of the unique values for the tables across both groups
        l_vec <- sapply(seen_q1_mat,length) + sapply(seen_q2_mat,length)
        l_tot <- sum(l_vec) # total number of tables across both groups

        # store the vector of weights W and build its mean value
        W0 <- c0/(c0+l_tot)
        W_vec <- l_vec/(c0+l_tot)
        mean_W0 <- mean_W0 + W0/R
        mean_W_vec <- mean_W_vec + W_vec/R

        # first component of the posterior variance of the baseline
        V0_1 <- V0_1 + (W0*int_diag_baseln + int_diag_tabs)/(c0 + l_tot + 1) - (1 - 1/(c0 + l_tot + 1))*(W0^2*int_cross_baseln + W0*int_cross_baseln_tabs + int_cross_tabs)

        # weight
        mean_W0 <- mean_W0 + W0/R
        mean_W_vec <- mean_W_vec + W_vec/R
        
        # compute the integral of the kernel k(x,y) w.r.t.
        # the empirical distribution of the dishes
        # according to the frequency of the tables across both groups both in x and in y
        int_cross_tabs <- quad.form(outer_mat_dishes,W_vec)

        # compute the integral of the kernel k(x,y) w.r.t.
        # the baseline measure in x,
        # the empirical distribution of the dishes
        # according to the frequency of the tables across both groups in y
        int_cross_baseln_tabs <- sum(colSums(outer_mat_baseln_dishes)*W_vec)/M

        # compute the integral of the kernel k(x,x) w.r.t
        # the empirical distribution of the dishes
        # according to the frequency of the tables across both groups in y
        int_diag_tabs <- sum(diag_vec_dishes*W_vec)
        
        if (n1_tot != 0) { # if at least one customer in group 1

          # compute the integral of the kernel k(x,y) w.r.t.
          # the empirical distribution of the dishes
          # according to the frequency of the tables across both groups in x,
          # the empirical distribution of the dishes
          # according to the frequency of the customers in group 1 in y
          int_cross_tabs_custs1 <- drop(quad.3form(outer_mat_dishes,W_vec,n1_vec))/n1_tot
        }
        else{ # if no customers in group 1

          # all the integrals involving the empirical distribution of the dishes
          # according to the frequency of the customers in group 1 are null
          int_cross_tabs_custs1 <- 0
        }
        
        # first and second components of the variance for group 1
        V1_1 <- V1_1 + (1 + 1/(c0 + l_tot + 1))/((c + n1 + 1)*(c + n1)^2)*((W0*int_diag_baseln + int_diag_tabs) - (W0^2*int_cross_baseln + W0*int_cross_baseln_tabs + int_cross_tabs))
        V1_2 <- V1_2 + c*n1/((c + n1 + 1)*(c + n1)^2)*((W0*int_diag_baseln + int_diag_tabs + int_diag_custs1) - (W0*int_cross_baseln_custs1 + int_cross_tabs_custs1))

        if (n2_tot != 0) { # if at least one customer in group 2

          # compute the integral of the kernel k(x,y) w.r.t.
          # the empirical distribution of the dishes
          # according to the frequency of the tables across both groups in x,
          # the empirical distribution of the dishes
          # according to the frequency of the customers in group 2 in y
          int_cross_tabs_custs2 <- drop(quad.3form(outer_mat_dishes,W_vec,n2_vec))/n2_tot
        }
        else{ # if no customers in group 2

          # all the integrals involving the empirical distribution of the dishes
          # according to the frequency of the customers in group 2 are null
          int_cross_tabs_custs2 <- 0
        }

        # first and second components of the variance for group 2
        V2_1 <- V2_1 + (1 + 1/(c0 + l_tot + 1))/((c + n2 + 1)*(c + n2)^2)*((W0*int_diag_baseln + int_diag_tabs) - (W0^2*int_cross_baseln + W0*int_cross_baseln_tabs + int_cross_tabs))
        V2_2 <- V2_1 + c*n2/((c + n2 + 1)*(c + n2)^2)*((W0*int_diag_baseln + int_diag_tabs + int_diag_custs2) - (W0*int_cross_baseln_custs2 + int_cross_tabs_custs2))
        
        # Gibbs update
        seen_q1_mat <- gibbs_tabs(l_vec,seen_q1_mat,c0=c0,c=c)
        l_vec <- sapply(seen_q1_mat,length) + sapply(seen_q2_mat,length)
        seen_q2_mat <- gibbs_tabs(l_vec,seen_q2_mat,c0=c0,c=c)
      }
    }
    else{
      return(hDP_corr_tabs(baseln = baseln,N_baseln = N_baseln,seen = matrix(numeric(0),ncol=4),c=c,c0=c0))
    }
    # compute the integral of the kernel k(x,y) w.r.t.
    # the empirical distribution of the dishes
    # according to the mean frequency of the tables across both groups both in x and in y
    int_cross_mean_tabs <- quad.form(outer_mat_dishes,mean_W_vec)

    # compute the integral of the kernel k(x,y) w.r.t.
    # the baseline measure in x,
    # the empirical distribution of the dishes
    # according to the mean frequency of the tables across both groups in y
    int_cross_baseln_mean_tabs <- sum(colSums(outer_mat_baseln_dishes)*mean_W_vec)/M

    # second component of the variance of the baseline
    V0_2 <- mean_W0^2*int_cross_baseln + mean_W0*int_cross_baseln_mean_tabs + int_cross_mean_tabs
    
    # variances of the two groups
    var1 <- V1_1 + V1_2 + V1_3 + c^2/(c + n1)^2*(V0_1 - V0_2)
    var2 <- V2_1 + V2_2 + V2_3 + c^2/(c + n2)^2*(V0_1 - V0_2)

    # covariance
    cov <- c^2/((c + n1)*(c + n2))*(V0_1 - V0_2)
    
    return(cov/sqrt(var1*var2))
}
