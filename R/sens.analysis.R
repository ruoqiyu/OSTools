sens.analysis<-function(y=NULL,D=NULL,Tobs=NULL,Gamma=1,tau=0,alternative='two-sided',method='m',
                        inner=0,trim=3,lambda=1/2,weight.par=c(1,1,1)){
  stopifnot((0 <= inner) & (inner <= trim))
  stopifnot((lambda > 0) & (lambda < 1))
  stopifnot(Gamma >= 1)
  m=weight.par[1]
  m1=weight.par[2]
  m2=weight.par[3]
  stopifnot((m1 >= 1) & (m2 >= m1) & (m >= m2))
  
  vc <- (sum(is.na(as.vector(y)))) > 0
  if (vc) warning("y cannot include NAs.")
  stopifnot(!vc)
  
  if (tau!=0) y[,1]<-y[,1]-tau
  
  if (method=='mcnemar'){
    if (is.null(D)){
      stopifnot(ncol(y)==2)
      D=sum(y[,1]!=y[,2])
    }
    if (is.null(Tobs)){
      stopifnot(ncol(y)==2)
      Tobs=sum(y[,1]>y[,2])
    }
    p.positive=Gamma/(1+Gamma);
    p.negative=1/(1+Gamma);
    lowerbound=1-pbinom(Tobs-1,D,p.negative);
    upperbound=1-pbinom(Tobs-1,D,p.positive);
    res=list(lowerbound=lowerbound,upperbound=upperbound);
    return(list(pval=res,statistic=Tobs))
  }
  
  if (method=='signedrank'){
    stopifnot(ncol(y)==2)
    dif=y[,1]-y[,2]
    rk=rank(abs(dif));
    s1=1*(dif>0);
    s2=1*(dif<0);
    W=sum(s1*rk);
    Eplus=sum((s1+s2)*rk*Gamma)/(1+Gamma);
    Eminus=sum((s1+s2)*rk)/(1+Gamma);
    V=sum((s1+s2)*rk*rk*Gamma)/((1+Gamma)^2);
    Dplus=(W-Eplus)/sqrt(V);
    Dminus=(W-Eminus)/sqrt(V);
    if (alternative=='greater') res=list(lowerbound=1-pnorm(Dminus),upperbound=1-pnorm(Dplus))
    else if (alternative=='less') res=list(lowerbound=pnorm(Dminus),upperbound=pnorm(Dplus))
    else if (alternative=='two-sided') res=list(lowerbound=2*pnorm(abs(Dminus)),upperbound=2*pnorm(abs(Dplus)))
    return(list(pval=res,statistic=W,deviate=list(Dminus=Dminus,Dplus=Dplus)))
  }
  
  if (method=='alignedrank'){
    outcome=as.vector(t(y))
    matchedset=rep(1:nrow(y),each=ncol(y))
    treated=rep(c(1,rep(0,ncol(y)-1)),nrow(y))
    # Compute means in each matched set
    matchedset.mean=tapply(outcome,matchedset,mean);
    # Compute residuals
    matchedset.mean.expand=matchedset.mean[matchedset];
    resids=outcome-matchedset.mean.expand;
    # Rank the residuals
    rankresids=rank(resids);
    # Test statistics = Sum of residuals in treatment group
    teststat=sum(rankresids[treated==1]);
    # Compute mu.i.max and sigma.i.max.sq in each matched set i
    # Assumes matched sets are labeled 1,...,I
    nomatchedset=length(unique(matchedset));
    mu.i.max=rep(0,nomatchedset);
    visq=rep(0,nomatchedset);
    for (i in 1:nomatchedset){
      ranks.matchedseti=rankresids[matchedset==i];
      notreated.matchedseti=sum(treated[matchedset==i]);
      if (notreated.matchedseti==1){
        sort.ranks.matchedseti=sort(ranks.matchedseti);
        ni=length(ranks.matchedseti);
        muia=rep(0,ni-1);
        viasq=rep(0,ni-1);
        for(j in 1:(ni-1)){
          muia[j]=(sum(sort.ranks.matchedseti[1:j])+Gamma*sum(sort.ranks.matchedseti[(j+1):ni]))/(j+Gamma*(ni-j));
          viasq[j]=(sum(sort.ranks.matchedseti[1:j]^2)+Gamma*sum(sort.ranks.matchedseti[(j+1):ni]^2))/(j+Gamma*(ni-j))-muia[j]^2;
        }
        mu.i.max[i]=max(muia);
        visq[i]=max(viasq[which(muia==max(muia))]);
      }
      if (notreated.matchedseti>1){
        sort.ranks.matchedseti=sort(ranks.matchedseti,decreasing=TRUE);
        ni=length(ranks.matchedseti);
        muia=rep(0,ni-1);
        viasq=rep(0,ni-1);
        totalranksum.matchedset=sum(ranks.matchedseti);
        for(j in 1:(ni-1)){
          muicontrol=(sum(sort.ranks.matchedseti[1:j])+Gamma*sum(sort.ranks.matchedseti[(j+1):ni]))/(j+Gamma*(ni-j));
          muia[j]=totalranksum.matchedset-muicontrol
          viasq[j]=(sum(sort.ranks.matchedseti[1:j]^2)+Gamma*sum(sort.ranks.matchedseti[(j+1):ni]^2))/(j+Gamma*(ni-j))-muicontrol^2;
        }
        mu.i.max[i]=max(muia);
        visq[i]=max(viasq[which(muia==max(muia))]);
      }
    }
    pval.u=1-pnorm((teststat-sum(mu.i.max))/sqrt(sum(visq)));
    
    # Compute mu.i.min and sigma.i.max.sq in each matched set i
    # Assumes matched sets are labeled 1,...,I
    nomatchedset=length(unique(matchedset));
    mu.i.min=rep(0,nomatchedset);
    visq=rep(0,nomatchedset);
    for (i in 1:nomatchedset){
      ranks.matchedseti=rankresids[matchedset==i];
      notreated.matchedseti=sum(treated[matchedset==i]);
      if (notreated.matchedseti==1){
        sort.ranks.matchedseti=sort(ranks.matchedseti);
        ni=length(ranks.matchedseti);
        muia=rep(0,ni-1);
        viasq=rep(0,ni-1);
        for (j in 1:(ni-1)){
          muia[j]=(sum(sort.ranks.matchedseti[1:j])+(1/Gamma)*sum(sort.ranks.matchedseti[(j+1):ni]))/(j+(1/Gamma)*(ni-j));
          viasq[j]=(sum(sort.ranks.matchedseti[1:j]^2)+(1/Gamma)*sum(sort.ranks.matchedseti[(j+1):ni]^2))/(j+(1/Gamma)*(ni-j))-muia[j]^2;
        }
        mu.i.min[i]=min(muia);
        visq[i]=max(viasq[which(muia==min(muia))]);
      }
      if (notreated.matchedseti>1){
        sort.ranks.matchedseti=sort(ranks.matchedseti,decreasing=TRUE);
        ni=length(ranks.matchedseti);
        muia=rep(0,ni-1);
        viasq=rep(0,ni-1);
        totalranksum.matchedset=sum(ranks.matchedseti);
        for (j in 1:(ni-1)){
          muicontrol=(sum(sort.ranks.matchedseti[1:j])+(1/Gamma)*sum(sort.ranks.matchedseti[(j+1):ni]))/(j+(1/Gamma)*(ni-j));
          muia[j]=totalranksum.matchedset-muicontrol
          viasq[j]=(sum(sort.ranks.matchedseti[1:j]^2)+(1/Gamma)*sum(sort.ranks.matchedseti[(j+1):ni]^2))/(j+(1/Gamma)*(ni-j))-muicontrol^2;
        }
        mu.i.min[i]=min(muia);
        visq[i]=max(viasq[which(muia==min(muia))]);
      }
    }
    pval.l=pnorm((teststat-sum(mu.i.min))/sqrt(sum(visq)));
    
    if (alternative=='greater') pval=pval.u
    else if (alternative=='less') pval=pval.l
    else if (alternative=='two-sided') pval=pval.u+pval.l
    return (list(pval=pval))
  }
  
  if (method=='m'){
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
    
    separable1k<-function(ymat,Gamma=1,alternative='two-sided'){
      stopifnot(0 == sum(is.na(as.vector(ymat))))
      n<-dim(ymat)[1]
      m<-dim(ymat)[2]
      o<-t(apply(ymat,1,sort))
      allmu<-matrix(NA,n,m-1)
      allsigma2<-matrix(NA,n,m-1)
      maxmu<-rep(-Inf,n)
      maxsig2<-rep(-Inf,n)
      for (j in 1:(m-1)){
        pr<-c(rep(1,j),rep(Gamma,m-j))/(j+((m-j)*Gamma))
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
      separable1k(newurks(ms,m=m,m1=m1,m2=m2),Gamma=Gamma,alternative=alternative)
    else if (m==1) 
      separable1k(ms,Gamma=Gamma,alternative=alternative)
    else warning('Consider choosing a diferent weight.par.')
  }
  
}
