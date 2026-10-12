sens.ci.full<-function(y,z,mset,Gamma=1,inner=0,trim=3,lambda=1/2,
                  alpha=0.05,alternative='two-sided'){
  alternative <- match.arg(alternative,c('two-sided','greater','less'))
  scalar <- function(x) is.numeric(x) && length(x)==1L && is.finite(x)
  if (!scalar(Gamma) || Gamma<1) stop("Gamma must be a finite scalar at least 1.")
  if (!scalar(alpha) || alpha<=0 || alpha>=1) stop("alpha must lie between 0 and 1.")
  if (!scalar(inner) || inner<0 || !is.numeric(trim) || length(trim)!=1L ||
      is.na(trim) || trim<inner || !scalar(lambda) || lambda<=0 || lambda>=1)
    stop("Score tuning requires 0 <= inner <= trim and 0 < lambda < 1.")
  if (!is.numeric(y) || !is.null(dim(y)) || !length(y)) stop("y must be a numeric response vector.")
  if (anyNA(z)) stop("z must contain nonmissing 0-1 treatment indicators.")

  stopifnot(length(y)==length(z), length(y)==length(mset),
            all(is.finite(y)), all(z %in% c(0,1)), !anyNA(mset))
  stratumsize=table(mset)
  mset01=table(z,mset)
  stopifnot(nrow(mset01)==2, all(mset01>0))
  stopifnot(all(apply(mset01,2,min)==1))
  nostratum=length(stratumsize)
  stratumlist=names(stratumsize)
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
  result <- sensitivityfull::senfmCI(y=ymat,treated1=treated1,gamma=Gamma,
      inner=inner,trim=trim,lambda=lambda,alpha=alpha,
      twosided=alternative=='two-sided',upper=alternative!='less')
  ci <- result$ConfidenceInterval
  if (!is.numeric(ci) || length(ci)!=2L || anyNA(ci) || ci[1]>ci[2])
    stop("The full-matching backend did not return an ordered confidence interval.")
  result
}
