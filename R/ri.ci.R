ri.ci=function(treated.r,control.r,alpha=0.05,alternative='greater',statistic='mean.diff',exact=FALSE,n.mc=1000,range=c(-100,100)){
  # Compute p-value-alpha/2 (reject for large T-tau-C)
  testtauci_lower=function(tau,treated.r,control.r,alpha,
                           statistic,exact,n.mc){
    ri.test(treated.r-tau,control.r,alternative='greater',statistic,exact,n.mc)$pval-alpha/(1+as.numeric(alternative=='two.sided'))
  }
  
  # Use bisection method to find lower endpoint of confidence interval
  lb=uniroot(testtauci_lower,range,
             treated.r,control.r,alpha,statistic,exact,n.mc)$root
  
  # Compute p-value-alpha/2 (reject for small T-tau-C)
  testtauci_upper=function(tau,treated.r,control.r,alpha,
                           statistic,exact,n.mc){ 
    ri.test(treated.r-tau,control.r,alternative='less',statistic,exact,n.mc)$pval-alpha/(1+as.numeric(alternative=='two.sided'))
  }
  
  # Use bisection method to find upper endpoint of confidence interval
  ub=uniroot(testtauci_upper,range,
             treated.r,control.r,alpha,statistic,exact,n.mc)$root
  if (alternative=='greater') ub=Inf
  if (alternative=='less') lb=-Inf
  list(lb=lb,ub=ub)
}

