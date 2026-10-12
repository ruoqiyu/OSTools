sens.ci<-function(y,Gamma=1,method='m',inner=0,trim=3,lambda=1/2,
                  weight.par=c(1,1,1),alpha=0.05,alternative='two-sided',
                  tol=NULL,interval=NULL){
  alternative <- match.arg(alternative, c('two-sided','greater','less'))
  method <- match.arg(method, c('m','signedrank','alignedrank'))
  scalar <- function(x) is.numeric(x) && length(x)==1L && is.finite(x)
  if ((!is.matrix(y) && !is.data.frame(y)) || nrow(y)<1L || ncol(y)<2L)
    stop("y must be a response matrix with at least one row and two columns.")
  y <- as.matrix(y)
  if (!is.numeric(y) || any(!is.finite(y))) stop("y must contain finite numeric responses.")
  if (method=='signedrank' && ncol(y)!=2L) stop("Signed ranks require matched pairs.")
  if (!scalar(Gamma) || Gamma<1) stop("Gamma must be a finite scalar at least 1.")
  if (!scalar(alpha) || alpha<=0 || alpha>=1) stop("alpha must lie between 0 and 1.")
  if (!scalar(inner) || inner<0 || !is.numeric(trim) || length(trim)!=1L ||
      is.na(trim) || trim<inner || !scalar(lambda) || lambda<=0 || lambda>=1)
    stop("Score tuning requires 0 <= inner <= trim and 0 < lambda < 1.")
  if (!is.numeric(weight.par) || length(weight.par)!=3L || any(!is.finite(weight.par)) ||
      any(weight.par!=floor(weight.par)) || weight.par[2]<1 ||
      weight.par[3]<weight.par[2] || weight.par[1]<weight.par[3])
    stop("weight.par must contain valid integer weights (m,m1,m2), with 1 <= m1 <= m2 <= m.")
  if (weight.par[3]!=weight.par[1]) stop("Redescending weights are not supported for confidence intervals.")
  if (alternative=='two-sided') alpha <- alpha/2
  if (is.null(interval)) interval <- c(min(y[,1])-max(y[,-1]),max(y[,1])-min(y[,-1]))
  if (!is.numeric(interval) || length(interval)!=2L || any(!is.finite(interval)) ||
      !is.finite(diff(interval)) || interval[2]<=interval[1])
    stop("interval must contain two finite, increasing endpoints.")
  if (is.null(tol)) tol <- diff(interval)/500000
  if (!scalar(tol) || tol<=0) stop("tol must be positive and finite.")
  if (method=='signedrank' && tol>diff(interval)) stop("Signed-rank grid spacing tol must not exceed the interval width.")
  root <- function(f, bounds, label, ...) {
    dots <- list(...)
    evaluate <- function(x) do.call(f,c(list(x),dots))
    values <- vapply(bounds,evaluate,numeric(1))
    if (any(!is.finite(values))) stop(paste("Cannot evaluate",label,"at the search endpoints; check responses and tuning."))
    if (all(values>0) || all(values<0))
      stop(paste("The search interval does not bracket",label,"; widen interval or check whether a finite bound is estimable."))
    stats::uniroot(evaluate,interval=bounds,f.lower=values[1],f.upper=values[2],
                   tol=tol,check.conv=TRUE)
  }

  if (method=='signedrank'){
    taugrid <- seq(interval[1],interval[2],tol)
    if (utils::tail(taugrid,1)<interval[2]) taugrid <- c(taugrid,interval[2])
    Dplus.taugrid=rep(0,length(taugrid));
    Dminus.taugrid=rep(0,length(taugrid));
    for(i in 1:length(taugrid)){
      tau=taugrid[i];
      deviates=sens.analysis(y,NULL,NULL,Gamma,tau,alternative,'signedrank',
                             inner,trim,lambda,weight.par)$deviate
      Dplus.taugrid[i]=deviates$Dplus;
      Dminus.taugrid[i]=deviates$Dminus;
    }
    
    critical <- stats::qnorm(alpha,lower.tail=FALSE)
    if (any(!is.finite(Dplus.taugrid)) || any(!is.finite(Dminus.taugrid)))
      stop("Signed-rank deviates are nonfinite; check responses and tuning.")
    need.lower <- alternative!='less'
    need.upper <- alternative!='greater'
    if (need.lower && Dplus.taugrid[1]<=critical)
      stop("The search interval truncates the lower confidence bound; widen interval or check whether a finite bound is estimable.")
    if (need.upper && utils::tail(Dminus.taugrid,1)>=-critical)
      stop("The search interval truncates the upper confidence bound; widen interval or check whether a finite bound is estimable.")
    accepted.lower <- Dplus.taugrid<=critical
    accepted.upper <- Dminus.taugrid>=-critical
    if ((need.lower && !any(accepted.lower)) || (need.upper && !any(accepted.upper)) ||
        (need.lower && need.upper && !any(accepted.lower & accepted.upper)))
      stop("No effects are retained on the supplied grid; widen interval or reduce tol.")
    lower.si <- if (!need.lower) -Inf else min(taugrid[accepted.lower])
    upper.si <- if (!need.upper) Inf else max(taugrid[accepted.upper])
    return(list(Confidence.Interval=c(lower.si,upper.si)))
  }
  
  if (method=='alignedrank'){
    cibound.func<-function(tau,Gamma,alternative="greater"){
      sens.analysis(y,Gamma=Gamma,tau=tau,alternative=alternative,method='alignedrank')$pval-alpha
    }
    
    # Confidence bounds
    lower <- if (alternative == 'less') -Inf else root(cibound.func,interval,'the lower confidence bound',Gamma=Gamma,alternative='greater')$root
    upper <- if (alternative == 'greater') Inf else root(cibound.func,interval,'the upper confidence bound',Gamma=Gamma,alternative='less')$root
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
    
    interval2<-c(-interval[2],-interval[1])
    vCI<-if (alternative == 'less') list(root=-Inf) else root(funcCI,interval,'the lower confidence bound',ymat=y)
    vEST<-root(funcEST,interval,'the minimum point estimate',ymat=y)
    vCI2<-if (alternative == 'greater') list(root=-Inf) else root(funcCI,interval2,'the upper confidence bound',ymat=-y)
    vEST2<-root(funcEST,interval2,'the maximum point estimate',ymat=-y)
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
