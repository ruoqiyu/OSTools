sens.ci.full<-function(y,z,mset,Gamma=1,inner=0,trim=3,lambda=1/2,
                  alpha=0.05,alternative='two-sided'){
  stratumsize=table(mset)
  mset01=table(z,mset)
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
  if (alternative=='two-sided') sensitivityfull::senfmCI(ymat,treated1,Gamma,inner,trim,lambda,alpha,T)   
  else if (alternative=='greater') sensitivityfull::senfmCI(ymat,treated1,Gamma,inner,trim,lambda,alpha,F,T) 
  else if (alternative=='less') sensitivityfull::senfmCI(ymat,treated1,Gamma,inner,trim,lambda,alpha,F,F) 
}
