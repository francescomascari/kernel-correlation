hDP_X2_to_mat <- function(
  seen_X2){
    X <- unique(seen_X2[,1])
    Ncusts1 <- sapply(X,function(val){sum(seen_X2[seen_X2[,2]==1,1] == val)})
    Ncusts2 <- sapply(X,function(val){sum(seen_X2[seen_X2[,2]==2,1] == val)})
    return(cbind(X=X,Ncusts1=Ncusts1,Ncusts2=Ncusts2))
}