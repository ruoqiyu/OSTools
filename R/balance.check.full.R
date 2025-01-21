balance.check.full=function(fdata,mdata,fz,mz,mset,fmissing=0,mmissing=0,parameter="ate"){
  stopifnot(dim(fdata)[2]==dim(mdata)[2])
  stopifnot(colnames(fdata)==colnames(mdata))
  
  msetcount=table(mset)
  mset01=table(mz,mset)
  nostratum=length(msetcount)
  stratumlist=names(msetcount)
  if(parameter=="ate"){
    weight=1/(.5/mset01[2,]+.5/mset01[1,])
  }
  if(parameter=="att"){
    weight=mset01[2,]
  }
  
  if (is.vector(fdata)){
    nvar=1
    if (is.null(nrow(fmissing))) fmissing=matrix(0,nrow=length(fdata),ncol=1)
    if (is.null(nrow(mmissing))) mmissing=matrix(0,nrow=length(mdata),ncol=1)
    fdata=as.matrix(fdata,ncol=1)
    mdata=as.matrix(mdata,ncol=1)
  }else{
    nvar=ncol(fdata)
    if (is.null(nrow(fmissing))) fmissing=matrix(0,nrow=nrow(fdata),ncol=ncol(fdata))
    if (is.null(nrow(mmissing))) mmissing=matrix(0,nrow=nrow(mdata),ncol=ncol(mdata))
  }
  
  smd.before=numeric(nvar)
  smd.after=numeric(nvar)
  for (k in 1:nvar){
    fx=fdata[,k]
    missingk=fmissing[,k]
    xtreated=fx[fz==1 & missingk==0];
    xcontrol=fx[fz==0 & missingk==0];
    var.xtreated=var(xtreated)
    var.xcontrol=var(xcontrol)
    combinedsd=sqrt(.5*(var.xtreated+var.xcontrol));
    smd.before[k]=(mean(xtreated)-mean(xcontrol))/combinedsd;
    
    x=mdata[,k]
    missing=mmissing[,k]
    diff.in.stratum=rep(0,nostratum);
    for(i in 1:nostratum){
      if(sum(mset==stratumlist[i] & mz==1 & missing==0)>0 & sum(mset==stratumlist[i] & mz==0 & missing==0)>0){
        diff.in.stratum[i]=mean(x[mset==stratumlist[i] & mz==1 & missing==0])-mean(x[mset==stratumlist[i] & mz==0 & missing==0]);
      }
      
    }
    smd.after[k]=(sum(weight*diff.in.stratum)/sum(weight))/combinedsd;
  }
  
  res=cbind(smd.before=smd.before,smd.after=smd.after)
  rownames(res)=colnames(fdata)
  res
}
