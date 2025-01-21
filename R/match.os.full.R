match.os.full<-function(z,dist,dat,ncontrol.min=0,ncontrol.max=Inf,exact=NULL,solver=''){
  #Check input
  stopifnot(is.data.frame(dat))
  stopifnot(is.matrix(dist))
  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))
  stopifnot(length(z)==(dim(dat)[1]))

  #sort input
  p=rep(1,length(z))
  if (is.null(exact)){
    o<-order(1-p)
  }else{
    o<-order(exact,1-p)
    exact<-exact[o]
  }
  
  z<-z[o]
  p<-p[o]
  dat<-dat[o,]
  
  #Must have treated first
  n=length(z)
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    o<-order(1-z)
    z<-z[o]
    dat<-dat[o,]
  }
  
  output=optmatch::fullmatch(dist,min.controls=ncontrol.min,max.controls=ncontrol.max,
                             data=dat,solver=solver)
  dat1=dat[!is.na(output),]
  output2=output[!is.na(output)]
  matchedset.index=substr(output2,start=3,stop=10)
  matchedset.index.numeric=as.numeric(matchedset.index)
  dat1$mset<-matchedset.index.numeric
  structure<-optmatch::stratumStructure(output)
  effsize<-optmatch::effectiveSampleSize(output)
  m<-list(data=dat1,structure=structure,effsize=effsize)
  m
}
