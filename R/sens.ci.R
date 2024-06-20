sens.ci<-function(y,gamma=1,method='m',inner=0,trim=3,lambda=1/2,
                  weight.par=c(1,1,1),alpha=0.05,alternative='two-sided',
                  tol=NULL,interval=NULL){
  m=weight.par[1]
  m1=weight.par[2]
  m2=weight.par[3]
  if (m2<m) warning("Redescending scores, m2<m, may not yield sensible confidence intervals and estimates")
  stopifnot(m2==m)
  
  stopifnot((alpha>0)&(alpha<1))
  if (alternative=='two-sided') alpha<-alpha/2
  if (method=='m'){
    funcCI<-function(tau,ymat){
      target<-alpha
      ntau<-length(tau)
      o<- rep(NA,ntau)
      for (i in 1:ntau){
        pp<-sens.analysis(ymat,gamma=gamma,tau=tau[i],alternative='greater',method=method,inner=inner,
                          trim=trim,lambda=lambda,weight.par=weight.par)$pval
        o[i]<-(pp-target)
      }
      o
    }
    
    funcEST<-function(tau,ymat){
      target<-0
      ntau<-length(tau)
      o<-rep(NA, ntau)
      for (i in 1:ntau){
        dev<-sens.analysis(ymat,gamma=gamma,tau=tau[i],alternative='greater',method=method,inner=inner,
                          trim=trim,lambda=lambda,weight.par=weight.par)$deviate
        o[i]<-(dev-target)
      }
      o
    }
    
    tr<-y[,1]
    ct<-as.vector(unlist(y[,-1]))
    mx<-max(tr)-min(ct)
    mn<-min(tr)-max(ct)
    if (is.null(interval)) interval<-c(mn,mx)
    else stopifnot(length(interval)==2)
    
    stopifnot(interval[2]>interval[1])
    interval2<-c(-interval[2],-interval[1])
    if (is.null(tol)) tol<-((max(interval)-min(interval)))/500000
    else stopifnot(tol>0)
    vCI<-stats::uniroot(funcCI,interval=interval,ymat=y,tol=tol)
    vEST<-stats::uniroot(funcEST,interval=interval,ymat=y,tol=tol)
    vCI2<-stats::uniroot(funcCI,interval=interval2,ymat=-y,tol=tol)
    vEST2<-stats::uniroot(funcEST,interval=interval2,ymat=-y,tol=tol)
    min.estimate<-vEST$root
    max.estimate<-(-vEST2$root)
    min.lowerCI<-vCI$root
    max.upperCI<-(-vCI2$root)
    PointEstimate<-c(min.estimate,max.estimate)
    names(PointEstimate)<-c("minimum","maximum")
    if (alternative=='greater') CI<-c(min.lowerCI,Inf)
    else if (alternative=='less') CI<-c(-Inf,max.upperCI)
    else if (alternative=='two-sided'){
      CI<-c(min.lowerCI,max.upperCI)
    }else{
      stop('alternative must be one of greater, less, and two-sided.')
    }
    names(CI)<-c("minimum","maximum")
    list(PointEstimate=PointEstimate,Confidence.Interval=CI)
  }else{
    warning('method needs to be \'m\'.')
  }
}
