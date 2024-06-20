sens.analysis<-function(y,gamma=1,tau=0,alternative='two-sided',method='m',
                        inner=0,trim=3,lambda=1/2,weight.par=c(1,1,1)){
  stopifnot((0 <= inner) & (inner <= trim))
  stopifnot((lambda > 0) & (lambda < 1))
  stopifnot(gamma >= 1)
  m=weight.par[1]
  m1=weight.par[2]
  m2=weight.par[3]
  stopifnot((m1 >= 1) & (m2 >= m1) & (m >= m2))
  
  vc <- (sum(is.na(as.vector(y)))) > 0
  if (vc) warning("y cannot include NAs.")
  stopifnot(!vc)
  
  if (tau!=0) y[,1]<-y[,1]-tau
  
  mscorev<-function(ymat,inner=0,trim=2.5,qu=0.5){
    ymat<-as.matrix(ymat)
    n<-dim(ymat)[1]
    m<-dim(ymat)[2]
    ou<-matrix(NA,n,m)
    one<-rep(1,m-1)
    difs<-array(NA,c(n,m,m-1))
    for (j in 1:m){
      difs[,j,]<-outer(as.vector(unlist(ymat[,j])),one,"*")-ymat[,-j]
    }
    ms<-as.vector(difs)
    if ((trim<Inf)|(inner>0)){
      hqu<-as.numeric(stats::quantile(abs(ms),qu,na.rm=TRUE))
      if (hqu>0){
        ms<-ms/hqu
        if ((trim<Inf)&(inner<trim)){
          ab<-pmin(1,pmax(0,(abs(ms)-inner))/(trim-inner))
        }else if ((trim<Inf)&(inner==trim)){
          ab<-1*(abs(ms)>inner)
        }else{
          ab<-pmax(0,abs(ms)-inner)
        }
        ms<-sign(ms)*ab
      }else{
        stop("Error: Scale factor is zero. Increase lambda.")
      }
    }
    ms<-array(ms,c(n,m,m-1))
    ms<-apply(ms,c(1,2),sum,na.rm=TRUE)
    ms[is.na(ymat)]<-NA
    colnames(ms)<-colnames(ymat)
    ni<-apply(!is.na(ymat),1,sum)
    use<-(ni>=2)&(!is.na(ms[, 1]))
    ms<-ms[use,]
    ni<-ni[use]
    
    ms<-ms/outer(ni,rep(1,m),"*")
    ms
  }
  
  ms<-mscorev(y,inner=inner,trim=trim,qu=lambda)
  
  if (method=='m'){
    newurks<-function(smat,m=1,m1=1,m2=1){
      rg<-apply(smat,1,max)-apply(smat,1,min)
      rk<-rank(rg)
      n<-length(rk)
      pk<-rk/n
      urk<-rep(0,n)
      for (l in m1:m2){
        urk<-urk+(l*choose(m,l)*(pk^(l-1))*((1-pk)^(m-l)))
      }
      for (j in 1:(dim(smat)[2])) {
        smat[,j]<-smat[,j]*urk
      }
      smat
    }
    
    separable1k<-function(ymat,gamma=1,alternative='two-sided'){
      stopifnot(0 == sum(is.na(as.vector(ymat))))
      n<-dim(ymat)[1]
      m<-dim(ymat)[2]
      o<-t(apply(ymat,1,sort))
      allmu<-matrix(NA,n,m-1)
      allsigma2<-matrix(NA,n,m-1)
      maxmu<-rep(-Inf,n)
      maxsig2<-rep(-Inf,n)
      for (j in 1:(m-1)){
        pr<-c(rep(1,j),rep(gamma,m-j))/(j+((m-j)*gamma))
        mu<-as.vector(o%*%pr)
        sigma2<-as.vector((o*o)%*%pr)-(mu*mu)
        chgmu<-(mu>maxmu)
        samemu<-(mu==maxmu)
        if (sum(chgmu)>0){
          maxmu[chgmu]<-mu[chgmu]
          maxsig2[chgmu]<-sigma2[chgmu]
        }
        if (sum(samemu)>0){
          maxsig2[samemu]<-pmax(sigma2[samemu],maxsig2[samemu])
        }
      }
      tstat<-as.vector(sum(ymat[,1]))
      expect<-sum(maxmu)
      vartotal<-sum(maxsig2)
      dev<-(tstat-expect)/sqrt(vartotal)
      if (alternative=='greater') pval<-1-stats::pnorm(dev)
      else if (alternative=='less') pval<-stats::pnorm(dev)
      else if (alternative=='two-sided') pval<-2*stats::pnorm(-abs(dev))
      else {
        pval=NULL
        warning ('alternative needs to be \'greater\', \'less\', or \'two-sided\'.')
      }
      list(pval=pval,deviate=dev,statistic=tstat, 
           expectation=expect,variance=vartotal)
    }
    
    if (m>1) 
      separable1k(newurks(ms,m=m,m1=m1,m2=m2),gamma=gamma,alternative=alternative)
    else if (m==1) 
      separable1k(ms,gamma=gamma,alternative=alternative)
    else warning('Consider choosing a different weight.par.')
  }else{
    warning('method needs to be \'m\'.')
  }
}
