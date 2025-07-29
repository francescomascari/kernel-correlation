my_seq_len <- function(n){
  if(n==0){
    return(c(0))
  }
  return(seq_len(n))
}