dilated.effect.matchedpair=function(treated.r,control.r,k,Delta0=0,alternative="greater",
                                    alpha=0.05,range=c(-100,100),
                                    grid=seq(0,10000,1),tol=0.0001){
  # Test of dilated treatment effect for matched pairs
  dilated.treateffect.matchedpair.test.func=function(Delta0,treated.r,control.r,k,
                                                     alternative="greater",
                                                     returntype="pval"){
    # Create vectors for Ri and Zi, and find total number in experiment and
    # number of treated.r subjects  
    Ri=c(treated.r,control.r);  
    Zi=c(rep(1,length(treated.r)),rep(0,length(control.r)));  
    N=length(Ri);  
    m=length(treated.r);    
    
    # Calculate adjusted responses and rho=r_{C(k)}  
    A=Ri-Zi*Delta0;  
    sorted.A=sort(A);  
    rho=sorted.A[k];    
    
    # q=1 if adjusted response>=rho, 0 otherwise  
    q=(A>=rho);
    qpaired=cbind(q[1:m],q[(m+1):(2*m)])  
    qsplus=apply(qpaired,1,sum)  
    no.discordant.pairs=sum(qsplus==1)    
    
    # Test statistic = # of discordant pairs in which treated.r has qs=1 
    teststat.obs= sum(qpaired[qsplus==1,1])    
    
    # For returning the p-value,  
    # p-value computed using McNemar’s test  
    if(returntype=="pval"& alternative=="less"){   
      returnval=pbinom(teststat.obs,no.discordant.pairs,.5)  
    }  
    if(alternative=="greater"){
      returnval=1-pbinom(teststat.obs-1,no.discordant.pairs,.5);  
    }    
    
    # For returning the test statistic minus its expected value 
    if(returntype=="teststat.minusev"){
      returnval=teststat.obs-no.discordant.pairs*.5;  
    }    
    returnval
  }
  
  # Search for endpoints of lower and upper bounds for confidence intervals
  pval.Delta0.func=function(Delta0,treated.r,control.r,k,alternative,alpha=0.05){
    dilated.treateffect.matchedpair.test.func(Delta0,treated.r,control.r,
                                              k,alternative)-alpha/2
  }
  
  dilated.ci=function(treated.r,control.r,k,range=c(-100,100),alpha=0.05){
    upper.ci.limit=uniroot(pval.Delta0.func,range,treated.r=treated.r,
                           control.r=control.r,k=k,alternative="less",alpha)$root;
    lower.ci.limit=uniroot(pval.Delta0.func,range,treated.r,
                           control.r,k=k,alternative="greater",alpha)$root;
    c(lower.ci.limit,upper.ci.limit)
  }
  
  # Hodges Lehmann estimate
  # Consider grid of values, find smallest value such that test statistic is 
  # less than its expectation, and largest value such that test statistic is 
  # greater than its expectation, and average these values
  dilated.hl.est=function(treated.r,control.r,k,grid=seq(0,10000,1),tol=0.0001){
    teststat.minus.ev.grid=rep(0,length(grid))  
    for(i in 1:length(teststat.minus.ev.grid)){
      teststat.minus.ev.grid[i]=
        dilated.treateffect.matchedpair.test.func(grid[i],treated.r,control.r,k,
                                                  returntype="teststat.minusev")  
    }  
    sup=max(grid[teststat.minus.ev.grid>tol]) 
    inf=min(grid[teststat.minus.ev.grid<(-tol)])  
    hl.est=(sup+inf)/2  
    hl.est
  }
  
  p=dilated.treateffect.matchedpair.test.func(Delta0,treated.r,control.r,k,alternative)
  ci=dilated.ci(treated.r,control.r,k,range,alpha)
  est=dilated.hl.est(treated.r,control.r,k,grid,tol)
  list(pvalue=p,conf.interval=ci,hl.estimate=est)
}
