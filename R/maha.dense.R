maha.dense<-function (z, X, exact=NULL, nearexact=NULL, penalty=100, matrix=FALSE){
  Xmatrix<-function(x){
    if (is.vector(x) || is.factor(x)) x<-matrix(x,nrow=length(z))

    if(is.data.frame(x) || is.character(x)){
      if(!is.data.frame(x)) x <- as.data.frame(x)
      X.chars <- which(plyr::laply(x, function(y) 'character' %in% class(y)))
      if(length(X.chars) > 0){
        for(i in X.chars){
          x[,i] <- factor(x[,i])

        }
      }
      #if some variables are factors convert to dummies
      X.factors <-  which(plyr::laply(x, function(y) 'factor' %in% class(y)))

      #handle missing data
      for(i in which(plyr::laply(x, function(y) any(is.na(y))))){
        if(i %in% X.factors){
          #for factors, make NA a new factor level
          x[,i] <- addNA(x[,i])
        }else{
          #for numeric/logical, impute means and add a new indicator for missingness
          x[[paste(colnames(x)[i],'NA', sep = '')]] <- is.na(x[,i])
          x[which(is.na(x[,i])),i] <- mean(x[,i], na.rm = TRUE)
        }
      }
      for(i in rev(X.factors)){
        dummyXi <- model.matrix(as.formula(
          paste('~',colnames(x)[i], '-1')),data=x)
        x <- cbind(x[,-i], dummyXi)
      }

    }else{
      #handle missing data
      for(i in c(1:ncol(x))){
        if(any(is.na(x[,i]))){
          x <- cbind(x,is.na(X[,i]))
          colnames(x)[ncol(x)] <- paste(colnames(X)[i],'NA', sep = '')
          x[which(is.na(x[,i])),i] <- mean(x[,i], na.rm = TRUE)
        }
      }

    }

    #get rid of columns that do not vary
    varying <- apply(x,2, function(y) length(unique(y)) > 1)
    x <- x[,which(varying),drop = FALSE]

    as.matrix(x)
  }

  X<-Xmatrix(X)
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


  for (j in 1:k) X[,j] <- rank(X[,j])
  cv <- cov(X)
  vuntied <- var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat) %*% cv %*% diag(rat)
#  LL<-chol(cv)
  icov <- MASS::ginv(cv)
  out <- matrix(NA, m, n-m)
  Xc <- X[z == 0,,drop=FALSE]
  Xt <- X[z == 1,,drop=FALSE]
  rownames(out) <- rownames(X)[z == 1]
  colnames(out) <- rownames(X)[z == 0]
#  for (i in 1:m) out[i, ] <- mvnfast::maha(Xc,t(as.matrix(Xt[i,])),LL,isChol=TRUE)
  for (i in 1:m) out[i, ] <- stats::mahalanobis(Xc,t(as.matrix(Xt[i,])),icov,inverted = T)
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
