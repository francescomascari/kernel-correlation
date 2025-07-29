log_un_prob_tabs2 <- function(
  n1_vec,
  n2_vec,
  l1_vec,
  l2_vec,
  c0,
  c){
    stopifnot(length(n1_vec) == length(n2_vec) && length(n1_vec)==length(l1_vec) && length(n2_vec)==length(l2_vec) )
    k <- length(n1_vec)
    n1_tot <- sum(n1_vec)
    n2_tot <- sum(n2_vec)
    
    l_vec <- l1_vec + l2_vec
    l_tot <- sum(l_vec)
    
    return(k*log(c0) + 2*lgamma(c)- lgamma(c+n1_tot) - lgamma(c+n2_tot) + l_tot*log(c) + lgamma(c0) - lgamma(c0+l_tot) + sum(log(abs(mapply(Stirling1,n1_vec,l1_vec)))) + sum(log(abs(mapply(Stirling1,n2_vec,l2_vec)))) + sum(lfactorial(l_vec - 1)))
}