sens.analysis.full<-function(y,z,mset,Gamma=1,method='m',alternative='two-sided',
                             tau=0,inner=0,trim=3,lambda=1/2){
  stratumsize=table(mset)
  mset01=table(z,mset)
  nostratum=length(stratumsize)
  stratumlist=names(stratumsize)
  if (method=='mh'){
    # Data for Mantel-Haenszel test
    data2x2vec=rep(0,nostratum)
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
