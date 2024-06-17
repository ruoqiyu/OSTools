maha.all<-function(X){
  Xmatrix<-function(x){
    if (is.vector(x) || is.factor(x)) x<-matrix(x,nrow=length(x))
    
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
  X<-as.matrix(X)
  n<-dim(X)[1]
  rownames(X)<-1:n
  k<-dim(X)[2]
  for (j in 1:k) X[,j]<-rank(X[,j])
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  icov<-MASS::ginv(cv)
  out<-matrix(NA,n,n)
  for (i in 1:n) out[i,]<-stats::mahalanobis(X,X[i,],icov,inverted=T)
  out
}
