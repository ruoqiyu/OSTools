sens.analysis.full<-function(y,z,mset,Gamma=1,method='m',alternative='two-sided',
                             tau=0,inner=0,trim=3,lambda=1/2){
  alternative <- match.arg(alternative, c('two-sided','greater','less'))
  method <- match.arg(method, c('m','mh'))
  if (!is.numeric(Gamma) || length(Gamma)!=1L || !is.finite(Gamma) || Gamma<1)
    stop("Gamma must be a finite scalar at least 1.")
  if (alternative == 'two-sided') {
    hi <- sens.analysis.full(y,z,mset,Gamma,method,'greater',tau,inner,trim,lambda)
    lo <- sens.analysis.full(y,z,mset,Gamma,method,'less',tau,inner,trim,lambda)
    result <- if (hi$pval <= lo$pval) hi else lo
    result$pval <- min(1,2*min(hi$pval,lo$pval))
    return(result)
  }
  stopifnot(length(y)==length(z), length(y)==length(mset),
            all(is.finite(y)), all(z %in% c(0,1)), !anyNA(mset))
  stratumsize=table(mset)
  mset01=table(z,mset)
  stopifnot(nrow(mset01)==2, all(mset01>0))
  if (method == 'm') stopifnot(all(apply(mset01,2,min)==1))
  nostratum=length(stratumsize)
  stratumlist=names(stratumsize)
  if (method=='mh'){
    stopifnot(all(y %in% c(0,1)), tau == 0)
    if (alternative == 'less') y <- 1-y
    # Data for Mantel-Haenszel test
    data2x2vec=rep(0,4*nostratum)
    index=0
    for(i in 1:nostratum){
      which.in.matched.set=which(mset==stratumlist[i])
      data2x2vec[index+1]=sum(z[which.in.matched.set]*y[which.in.matched.set])
      data2x2vec[index+2]=sum((1-z[which.in.matched.set])*y[which.in.matched.set])
      data2x2vec[index+3]=sum(z[which.in.matched.set]*(1-y[which.in.matched.set]))
      data2x2vec[index+4]=sum((1-z[which.in.matched.set])*(1-y[which.in.matched.set]))
      index=index+4
    }
    data.array=array(data2x2vec,c(2,2,nostratum))
    # Sensitivity analysis 
    # Sum the upper tail as a lower tail of failures, using positive-only
    # truncated convolution. FFT convolution can swamp very small tails.
    m1 <- data.array[1,1,]+data.array[1,2,]
    m2 <- data.array[2,1,]+data.array[2,2,]
    n <- data.array[1,1,]+data.array[2,1,]
    minimum <- pmax(0,n-m2)
    maximum <- pmin(n,m1)
    observed <- sum(data.array[1,1,])
    target <- sum(maximum)-observed
    probabilities <- 1
    for (i in seq_len(nostratum)) {
      x <- seq.int(minimum[i],maximum[i])
      logweights <- lchoose(m1[i],x)+lchoose(m2[i],n[i]-x)+
        (x-maximum[i])*log(Gamma)
      weights <- exp(logweights-max(logweights))
      failures <- rev(weights/sum(weights))
      size <- min(target+1,length(probabilities)+length(failures)-1)
      updated <- numeric(size)
      for (j in seq_len(min(length(failures),size))) {
        use <- seq_len(min(length(probabilities),size-j+1))
        updated[use+j-1] <- updated[use+j-1]+failures[j]*probabilities[use]
      }
      probabilities <- updated
    }
    pval <- min(1,sum(probabilities))
    list(pval=pval,data.array=data.array)
  }else if (method=='m'){
    treated1=rep(0,nostratum)
    ymat=matrix(rep(NA,nostratum*max(stratumsize)),nrow=nostratum)
    for(i in 1:nostratum){
      no.treated.in.stratum=mset01[2,i]
      no.control.in.stratum=mset01[1,i]
      treated.in.stratum=which(mset==stratumlist[i] & z==1)
      control.in.stratum=which(mset==stratumlist[i] & z==0)  
      if(no.treated.in.stratum==1){
        ymat[i,1]=y[treated.in.stratum]
        ymat[i,2:(no.control.in.stratum+1)]=y[control.in.stratum]
        treated1[i]=1
      }
      if(no.treated.in.stratum>1){
        ymat[i,1]=y[control.in.stratum]
        ymat[i,2:(no.treated.in.stratum+1)]=y[treated.in.stratum]
        treated1[i]=0
      }
    }
  
    treated1=as.logical(treated1)
    result <- sensitivityfull::senfm(ymat,treated1,Gamma,inner,trim,lambda,tau,alternative)
    # The backend deviate already represents the requested direction.
    result$pval <- stats::pnorm(result$deviate,lower.tail=FALSE)
    result   
  }
}
