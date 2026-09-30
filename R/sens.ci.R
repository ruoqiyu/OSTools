sens.ci<-function(y,Gamma=1,method='m',inner=0,trim=3,lambda=1/2,
                  weight.par=c(1,1,1),alpha=0.05,alternative='two-sided',
                  tol=NULL,interval=NULL){
  alternative <- match.arg(alternative, c('two-sided','greater','less'))
  method <- match.arg(method, c('m','signedrank','alignedrank'))
  m=weight.par[1]
  m1=weight.par[2]
  m2=weight.par[3]
  if (m2<m) warning("Redescending scores, m2<m, may not yield sensible confidence intervals and estimates")
  stopifnot(m2==m)
  
  stopifnot((alpha>0)&(alpha<1))
  if (alternative=='two-sided') alpha<-alpha/2
  
  if (method=='signedrank'){
    tr<-y[,1]
    ct<-as.vector(unlist(y[,-1]))
    mx<-max(tr)-min(ct)
    mn<-min(tr)-max(ct)
    if (is.null(interval)) interval<-c(mn,mx)
    else stopifnot(length(interval)==2)
    stopifnot(interval[2]>interval[1])
#    interval2<-c(-interval[2],-interval[1])
    if (is.null(tol)) tol<-((max(interval)-min(interval)))/500000
    else stopifnot(tol>0)
    
    taugrid=seq(interval[1],interval[2],tol);
    Dplus.taugrid=rep(0,length(taugrid));
    Dminus.taugrid=rep(0,length(taugrid));
    for(i in 1:length(taugrid)){
      tau=taugrid[i];
      deviates=sens.analysis(y,NULL,NULL,Gamma,tau,alternative,'signedrank',
                             inner,trim,lambda,weight.par)$deviate
      Dplus.taugrid[i]=deviates$Dplus;
      Dminus.taugrid[i]=deviates$Dminus;
    }
    
    critical <- stats::qnorm(1-alpha)
    lower.si <- if (alternative == 'less') -Inf else min(taugrid[Dplus.taugrid<=critical])
    upper.si <- if (alternative == 'greater') Inf else max(taugrid[Dminus.taugrid>=-critical])
    return(list(Confidence.Interval=c(lower.si,upper.si)))
  }
  
  if (method=='alignedrank'){
    cibound.func<-function(tau,Gamma,alternative="greater"){
      sens.analysis(y,Gamma=Gamma,tau=tau,alternative=alternative,method='alignedrank')$pval-alpha
    }
    
    # Confidence bounds
    if (is.null(interval)) interval <- c(min(y[,1])-max(y[,-1]), max(y[,1])-min(y[,-1]))
    if (is.null(tol)) tol <- diff(interval)/500000
    lower <- if (alternative == 'less') -Inf else uniroot(cibound.func,interval,Gamma=Gamma,alternative='greater',tol=tol)$root
    upper <- if (alternative == 'greater') Inf else uniroot(cibound.func,interval,Gamma=Gamma,alternative='less',tol=tol)$root
    return(list(Confidence.Interval=c(lower,upper)))
  }
  
  if (method=='m'){
    funcCI<-function(tau,ymat){
      target<-alpha
      ntau<-length(tau)
      o<- rep(NA,ntau)
      for (i in 1:ntau){
        pp<-sens.analysis(ymat,Gamma=Gamma,tau=tau[i],alternative='greater',method=method,inner=inner,
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
        dev<-sens.analysis(ymat,Gamma=Gamma,tau=tau[i],alternative='greater',method=method,inner=inner,
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
    vCI<-if (alternative == 'less') list(root=-Inf) else stats::uniroot(funcCI,interval=interval,ymat=y,tol=tol)
    vEST<-stats::uniroot(funcEST,interval=interval,ymat=y,tol=tol)
    vCI2<-if (alternative == 'greater') list(root=-Inf) else stats::uniroot(funcCI,interval=interval2,ymat=-y,tol=tol)
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
  }
}
