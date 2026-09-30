sens.ci.iv<-function(encouraged.treatment,unencouraged.treatment,
                     encouraged.response,unencouraged.response,
                     Gamma=1,alpha=0.05,tol=0.01,interval=c(-10,10)){
  
  stopifnot((alpha>0)&(alpha<1))
  
  sens.interval.iv.deviates=function(beta,encouraged.treatment,unencouraged.treatment,
                                     encouraged.response,unencouraged.response,Gamma){
    adjusted.response.encouraged=encouraged.response-beta*encouraged.treatment;
    adjusted.response.unencouraged=unencouraged.response-beta*unencouraged.treatment;
    diff=adjusted.response.encouraged-adjusted.response.unencouraged;
    rk=rank(abs(diff));
    s1=1*(diff>0);
    s2=1*(diff<0);
    W=sum(s1*rk);
    Eplus=sum((s1+s2)*rk*Gamma)/(1+Gamma);
    Eminus=sum((s1+s2)*rk)/(1+Gamma);
    V=sum((s1+s2)*rk*rk*Gamma)/((1+Gamma)^2);
    Dplus=(W-Eplus)/sqrt(V);
    Dminus=(W-Eminus)/sqrt(V);
    list(Dplus=Dplus,Dminus=Dminus);
  }
  
  betagrid=seq(interval[1],interval[2],tol)
  Dplus.betagrid=rep(0,length(betagrid));
  Dminus.betagrid=rep(0,length(betagrid));
  for(i in 1:length(betagrid)){
    beta=betagrid[i];
    deviates=sens.interval.iv.deviates(beta,encouraged.treatment,unencouraged.treatment,
                                       encouraged.response,unencouraged.response,Gamma);
    Dplus.betagrid[i]=deviates$Dplus;
    Dminus.betagrid[i]=deviates$Dminus;
  }
  
  lower.si=min(betagrid[Dplus.betagrid<=qnorm(1-alpha/2)]);
  upper.si=max(betagrid[Dminus.betagrid>=qnorm(alpha/2)]);
  
  c(lower.si,upper.si)
}
