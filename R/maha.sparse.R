maha.sparse<-function(z,X,p=rep(1,length(z)),caliper=1,stdev=FALSE,constant=NULL,ncontrol=1,
                      exact=NULL,nearexact=NULL,penalty=100,subX=NULL,ties.all=TRUE){


  #Check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  stopifnot(all((z==1)|(z==0)))
  nobs<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(length(z)==length(p))
  X<-.os_covariates(X)

  if (is.factor(nearexact)){
    levels(nearexact)=1:nlevels(nearexact)
    nearexact<-as.integer(nearexact)
  }

  stopifnot(length(z)==(dim(X)[1]))

  if (!is.null(exact)){
    stopifnot(length(exact)==length(z))
    tb<-table(z,exact)
    if (!all(ncontrol*tb[2,]<=tb[1,])){
      if (is.null(subX)|| any(exact!=subX)){
        warning("An exact match for exact is infeasible for every caliper.")
        return(NULL)
      }
    }
  }

  if (length(caliper)==1) stopifnot(caliper>=0)

  # Standardize p using an equally weighted pooled variance
  if (stdev){
    v1<-stats::var(p[z==1])
    v2<-stats::var(p[z==0])
    sp<-sqrt((v1+v2)/2)
    caliper<-caliper*sp
  }

  n<-dim(X)[1]

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
  if (!is.null(subX)) subX<-subX[o]
  if (is.vector(X)) X<-matrix(X,length(X),1)

  #Must have treated first
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    o<-order(1-z)
    z<-z[o]
    p<-p[o]
    X<-X[o,,drop=FALSE]
    if (!is.null(exact)) exact<-exact[o]
    if (!is.null(nearexact)) nearexact<-nearexact[o]
  }

  if (is.vector(X)) X<-matrix(X,ncol=1)

  ids<-1:n
  k<-dim(X)[2]
  m<-sum(z)

  for (j in seq_len(k)) X[, j]<-rank(X[, j])
  cv <- if (k == 0L) matrix(0,0,0) else stats::cov(X)
  rat <- sqrt(stats::var(seq_len(n))/diag(cv))
  cv <- diag(rat,nrow=length(rat)) %*% cv %*% diag(rat,nrow=length(rat))
  icov <- if (k == 0L) cv else MASS::ginv(cv)
  X0 <- X[z==0,,drop=FALSE]
  X1 <- X[z==1,,drop=FALSE]
  cids <-ids[z==0]
  p0<-p[z==0]
  p1<-p[z==1]
  if (!is.null(exact)){
    exact0<-exact[z==0]
    exact1<-exact[z==1]
  }
  if (!is.null(nearexact)){
    nearex0<-nearexact[z==0]
    nearex1<-nearexact[z==1]
  }

  if (length(caliper)==1) caliper<-c(-caliper,caliper)
  if (length(caliper)==2) caliper<-sort(caliper)

  edgenum<-function(z,p,caliper,constant=NULL,exact=NULL,ties.all=TRUE){

    stopifnot(is.vector(z))
    stopifnot(all((z==1)|(z==0)))
    stopifnot(caliper>=0)
    if (!is.null(constant)){
      if (!is.integer(constant)) constant<-floor(constant)
    }

    if (!is.null(exact)){
      stopifnot(length(exact)==length(z))
    }

    n<-length(z)

    #Must have treated first
    if(!(min(z[1:(n-1)]-z[2:n])>=0)){
      o<-order(1-z)
      z<-z[o]
      p<-p[o]
      if (!is.null(exact)) exact<-exact[o]
    }

    ids<-1:n
    m<-sum(z)

    cids <-ids[z==0]
    p0<-p[z==0]
    p1<-p[z==1]
    if (!is.null(exact)){
      exact0<-exact[z==0]
      exact1<-exact[z==1]
    }


    result=0

    for (i in 1:m){
      #use caliper
      d<-abs(p1[i]-p0)
      who<-d<=caliper
      #use exact
      if (!is.null(exact)) who<-who&(exact1[i]==exact0)
      num<-sum(who)
      #use constant
      if (!is.null(constant)){
        if(num>constant){
          #        who<-which(who)[order(d[who])[1:constant]]
          #        num<- constant
          whoi<-which(who)[rank(d[who],ties.method='min')<=constant]
          if (ties.all || length(whoi)==constant){
            who<-whoi
            num<-length(who)
          }else{
            num<- constant
          }
        }
      }
      result=result+num
    }
    result
  }

  edgen<-edgenum(z,p,max(abs(caliper)),constant,exact,ties.all)
  distance<-numeric(edgen)
  start_node<-numeric(edgen)
  end_node<-numeric(edgen)
  if (!is.null(nearexact)){
    nearex<-numeric(edgen)
  }else{
    nearex<-c()
  }

  left<-numeric(m)
  right<-numeric(m)

  current<-0

  for (i in 1:m){
    #use caliper
    d<-p1[i]-p0
    who<-(d<=caliper[2])&(d>=caliper[1])
    #use exact
    d<-abs(p1[i]-p0)
    if (!is.null(exact)) who<-who&(exact1[i]==exact0)
    if (!any(who)) return(NULL)
    num<-sum(who)
    #use constant
    nexti<-ncontrol
    if (!is.null(constant)){
      if(num>constant){
        #        who<-which(who)[order(d[who])[1:constant]]
        #        num<- constant
        whoi<-which(who)[rank(d[who],ties.method='min')<=constant]
        if (ties.all || length(whoi)==constant){
          who<-whoi
          num<-length(who)
        }else if (constant==ncontrol){
          whoi<-which(who)[order(d[who])[1:ncontrol]]
          if ((cids[whoi][1] %in% left) && (cids[whoi][constant] %in% right)){
            nexti<-nexti+1
            if (nexti+constant-1>n-m) return(NULL)
            else{
              nextclose<-all(who[order(d[who])[1:constant]]==who[order(d[who])[nexti:(nexti+constant-1)]])
              whoi<-which(who)[order(d[who])[nexti:(nexti+constant-1)]]
              while(nextclose && (cids[whoi][1] %in% left) && (cids[whoi][constant] %in% right) && (nexti+constant-1<=n-m)){
                nexti<-nexti+1
                if (nexti+constant-1>n-m) return(NULL)
                whoi<-which(who)[order(d[who])[nexti:(nexti+constant-1)]]
              }
            }
          }
          who<-whoi
          num<- constant
          left[i]<-min(cids[who])
          right[i]<-max(cids[who])
        }else if(any(d[whoi]!=max(d[whoi]))){
          closea<-whoi[which(d[whoi]!=max(d[whoi]))]
          closeb<-which(who)[which(d[who]==max(d[whoi]))]
          closeb<-closeb[order(pmin(abs(closeb-min(closea)),abs(closeb-max(closea))))]
          closeb<-closeb[1:(constant-length(closea))]
          who<-c(closea,closeb)
          num<- constant
        }else{
          who<-which(who)[1:constant]
          num<- constant
        }
      }
    }

    cc<-X0[who,,drop=FALSE]
    if (num==1) cc<-matrix(cc,nrow=1)
    else if (is.vector(cc)) cc<-matrix(cc,nrow=num)
    tt<-t(as.matrix(X1[i,]))
    distancei<-if (k == 0L) rep(0,num) else stats::mahalanobis(cc,tt,icov,inverted=TRUE)

    if (num>0){
      distance[(current+1):(current+num)]<-distancei
      start_node[(current+1):(current+num)]<-rep(i, num)
      end_node[(current+1):(current+num)]<-cids[who]

      #use nearexact
      if (!is.null(nearexact)) nearex[(current+1):(current+num)]<-nearex1[i]!=nearex0[who]
    }

    current<-current+num
  }
  if (!is.null(nearexact)){
    out<-list(d=(distance+nearex*penalty)[seq_len(current)],start=start_node[seq_len(current)],end=end_node[seq_len(current)])
  }else{
    out<-list(d=distance[seq_len(current)],start=start_node[seq_len(current)],end=end_node[seq_len(current)])
  }
  out
}
