sens.analysis.full<-function(y,z,mset,Gamma=1,method='m',alternative='two-sided',
                             tau=0,inner=0,trim=3,lambda=1/2){
  alternative <- match.arg(alternative, c('two-sided','greater','less'))
  method <- match.arg(method, c('m','mh'))
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
    pval=sensitivity2x2xk::mh(data.array,Gamma=Gamma)$pval
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
    sensitivityfull::senfm(ymat,treated1,Gamma,inner,trim,lambda,tau,alternative)   
  }
}
