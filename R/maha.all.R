maha.all<-function(X){

  
  X<-.os_covariates(X)
  X<-as.matrix(X)
  n<-dim(X)[1]
  rownames(X)<-1:n
  k<-dim(X)[2]
  for (j in seq_len(k)) X[,j]<-rank(X[,j])
  if (k == 0L) return(matrix(0,n,n))
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat,nrow=length(rat))%*%cv%*%diag(rat,nrow=length(rat))
  icov<-MASS::ginv(cv)
  out<-matrix(NA,n,n)
  for (i in 1:n) out[i,]<-stats::mahalanobis(X,X[i,],icov,inverted=T)
  out
}
