hDP_XT2_to_mat <- function(
  seen_XT2){
    X <- unique(seen_XT2[,1])
    Ntabs <- sapply(X,function(val){length(unique(seen_XT2[seen_XT2[,1] == val,2]))})
    Ncusts1 <- sapply(X,function(val){sum(seen_XT2[seen_XT2[,3]==1,1] == val)})
    Ncusts2 <- sapply(X,function(val){sum(seen_XT2[seen_XT2[,3]==2,1] == val)})
    return(data.frame(X=X,Ntabs=Ntabs,Ncusts1=Ncusts1,Ncusts2=Ncusts2))
}