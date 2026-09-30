maha.dense<-function (z, X, exact=NULL, nearexact=NULL, penalty=100, matrix=FALSE){


  stopifnot(is.vector(z), all(z %in% c(0,1)), all(c(0,1) %in% z))
  X<-.os_covariates(X)
  stopifnot(nrow(X) == length(z))
  n <- dim(X)[1]
  rownames(X) <- 1:n
  k <- dim(X)[2]
  m <- sum(z)

  p=rep(1,length(z))
  #sort input
  if (is.null(exact)){
    o<-order(1-p)
  }else{
    o<-order(exact,1-p)
    exact<-exact[o]
  }

  z<-z[o]
  p<-p[o]
  X<-X[o,,drop=FALSE]
  if (!is.null(nearexact)) nearexact<-nearexact[o]

  #Must have treated first
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    o<-order(1-z)
    z<-z[o]
    p<-p[o]
    X<-X[o,,drop=FALSE]
    if (!is.null(exact)) exact<-exact[o]
    if (!is.null(nearexact)) nearexact<-nearexact[o]
  }


  for (j in seq_len(k)) X[,j] <- rank(X[,j])
  cv <- if (k == 0L) matrix(0,0,0) else cov(X)
  vuntied <- var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat,nrow=length(rat)) %*% cv %*% diag(rat,nrow=length(rat))
#  LL<-chol(cv)
  icov <- if (k == 0L) cv else MASS::ginv(cv)
  out <- matrix(NA, m, n-m)
  Xc <- X[z == 0,,drop=FALSE]
  Xt <- X[z == 1,,drop=FALSE]
  rownames(out) <- rownames(X)[z == 1]
  colnames(out) <- rownames(X)[z == 0]
#  for (i in 1:m) out[i, ] <- mvnfast::maha(Xc,t(as.matrix(Xt[i,])),LL,isChol=TRUE)
  for (i in seq_len(m)) out[i, ] <- if (k == 0L) 0 else stats::mahalanobis(Xc,t(as.matrix(Xt[i,])),icov,inverted = TRUE)
  if (!is.null(exact)){
    dif <- outer(exact[z == 1], exact[z == 0], "!=")
    out[dif] <- Inf
  }

  if (!is.null(nearexact)){
    dif <- outer(nearexact[z == 1], nearexact[z == 0], "!=")
    out <- out + dif * penalty
  }

  if (matrix) return (out)
  else{
    distance<-t(out)
    dim(distance)<-c(1,m*(n-m))
    distance<-as.vector(distance)
    start<-rep(1:m,each=n-m)
    end<-rep((m+1):n,m)
    d0<-distance
    distance<-distance[which(d0<Inf)]
    start<-start[which(d0<Inf)]
    end<-end[which(d0<Inf)]
    return (list(d=distance,start=start,end=end))
  }
}
