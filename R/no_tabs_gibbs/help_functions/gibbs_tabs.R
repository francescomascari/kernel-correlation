gibbs_tabs <- function(
  l_vec,
  seen_q_vec,
  c0 = 1,
  c = 1){
    n_dishes <- length(l_vec)
    l_tot <- sum(l_vec)
    
    for(n_indx in seq_len(n_dishes)){
      start_q_vec <- seen_q_vec[n_indx][[1]]
      q_vec <- start_q_vec
      
      len_q <- length(start_q_vec)
      
      for(q_indx in seq_len(len_q)){
        q <- start_q_vec[q_indx]
        for(i in seq_len(q)){
          
          if(q_vec[q_indx] == 1){
            l_vec[n_indx] <- l_vec[n_indx] - 1
            l_tot <- l_tot - 1
          }
          q_vec[q_indx] <- q_vec[q_indx] - 1
          
          probs <- update_q_probs(q_vec,l_tot,l_vec[n_indx],c0=c0,c=c)
          new_indx <- sample(0:length(q_vec), size = 1,prob = probs)
          
          if(new_indx == 0){
            q_vec <- c(q_vec,1)
            l_vec[n_indx] <- l_vec[n_indx] + 1
            l_tot <- l_tot + 1
          }
          else{
            q_vec[new_indx] <- q_vec[new_indx] + 1
          }
        }
      }
      seen_q_vec[n_indx][[1]] <- q_vec[q_vec != 0]
    }
    return(seen_q_vec)
}