update_q_probs <- function(
  q_vec_k,
  l_tot,
  l_k,
  c=1,
  c0=1){
    len_q <- length(q_vec_k)
    
    if(l_k == 0){
      return(c(1,rep(0,len_q)))
    }
    
    un_probs <- rep(numeric(0),len_q + 1)
    un_probs[1] <- log(c) - log(c0 + l_tot) + log(l_k)
    for(i in seq_len(len_q)){
      un_probs[i+1] <- log(q_vec_k[i]) #2*log(q_vec_k[i])
    }
    
    probs <- exp(un_probs - logSumExp(un_probs))
    
    return(probs)
}