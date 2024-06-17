balance.check<-function(fdata=NULL,mdata=NULL,fz=NULL,mz=NULL,measure='smd',s.dist=1000){
  if (is.null(fdata) & is.null(mdata))
    stop('At least one of fdata and mdata needs to be entered.')

  if ((!is.null(fdata)) & (!is.null(mdata))){
    stopifnot(dim(fdata)[2]==dim(mdata)[2])
    stopifnot(colnames(fdata)==colnames(mdata))
    stopifnot(all(fz%in%c(0,1)))
    stopifnot(all(mz%in%c(0,1)))

    if (is.vector(fdata)) fdata<-matrix(fdata,ncol=1)
    if (is.vector(mdata)) mdata<-matrix(mdata,ncol=1)
    fd0<-fdata[which(fz==0),]
    fd1<-fdata[which(fz==1),]
    md0<-mdata[which(mz==0),]
    md1<-mdata[which(mz==1),]
    if (is.vector(fd0)) fd0<-matrix(fd0,length(fd0),1)
    if (is.vector(md0)) md0<-matrix(md0,length(md0),1)
    if (is.vector(fd1)) fd1<-matrix(fd1,length(fd1),1)
    if (is.vector(md1)) md1<-matrix(md1,length(md1),1)

    if (measure=='smd'){
      varf0<-apply(fd0,2,var,na.rm=TRUE)
      varf1<-apply(fd1,2,var,na.rm=TRUE)
      meanf0<-apply(fd0,2,mean,na.rm=TRUE)
      meanf1<-apply(fd1,2,mean,na.rm=TRUE)
      meanm0<-apply(md0,2,mean,na.rm=TRUE)
      meanm1<-apply(md1,2,mean,na.rm=TRUE)
      smdf<-(meanf1-meanf0)/sqrt((varf0+varf1)/2)
      smdm<-(meanm1-meanm0)/sqrt((varf0+varf1)/2)

      r<-cbind(meanf1,meanm0,meanf0,smdm,smdf)
      rownames(r)<-colnames(fdata)
      colnames(r)<-c('Treated Mean','Control Match Mean','Control All Mean',
                     'After SMD','Before SMD')
      r
    }else if (measure=='t.test'){
      nvar=ncol(fdata)
      meanf0<-apply(fd0,2,mean,na.rm=TRUE)
      meanf1<-apply(fd1,2,mean,na.rm=TRUE)
      meanm0<-apply(md0,2,mean,na.rm=TRUE)
      meanm1<-apply(md1,2,mean,na.rm=TRUE)
      t.test.pval.before=rep(0,nvar)
      t.test.pval.after=rep(0,nvar)
      for (i in 1:nvar){
        t.test.pval.before[i]=
          stats::t.test(fd1[,i],fd0[,i])$p.value
        t.test.pval.after[i]=
          stats::t.test(md1[,i],md0[,i])$p.value
      }
      r=cbind(meanf1,meanm0,meanf0,t.test.pval.after,t.test.pval.before)
      rownames(r)<-colnames(fdata)
      colnames(r)<-c('Treated Mean','Control Match Mean','Control All Mean',
                     'After P-value','Before P-value')
      r
    }else if (measure=='crossmatch'){
      D=maha.all(mdata)
      plainmatrix<-as.matrix(s.dist*max(1/min(D[D>0]),1)*D)
      diag(plainmatrix) <- 0
      mdm<-nbpMatching::distancematrix(plainmatrix)
      res<-nbpMatching::nonbimatch(mdm)

      mt_min<-pmin(as.numeric(res$matches$Group1.Row),
                   as.numeric(res$matches$Group2.Row))
      mt_max<-pmax(as.numeric(res$matches$Group1.Row),
                   as.numeric(res$matches$Group2.Row))
      z0<-mz[(mt_min>0)&(mt_max<=length(mz))]
      mt0<-factor(mt_min[(mt_min>0)&(mt_max<=length(mz))])
      tab<-table(factor(z0),mt0)
      a1<-sum(tab[1,]==1)

      crossmatchdist<-function(N,n){
        if (bigN%%2 == 1){
          stop("The number of subjects, bigN, should be even")
          return(NA)
        }
        I<-N/2
        cmd<-c()
        for (a1 in 0:I){
          a2 <- (n-a1)/2
          if ((floor(a2)==a2)&(a2>=0)){
            a0<-I-(a1+a2)
            if (a0>=0){
              pr<-factorial(I)/choose(N,n)
              pr<-pr*(2^a1)/(factorial(a0)*factorial(a1)*factorial(a2))
              cmd<-cbind(cmd,c(a0,a1,a2,pr))
            }
          }
        }
        cmd <- rbind(cmd,cumsum(cmd[4,]))
        cmd
      }

      bigN<-length(z0)
      smalln<-sum(z0)
      if (bigN<340){
        dist<-crossmatchdist(bigN,smalln)
        pval<-dist[5,dist[2,]==a1]
      }else{
        pval<-NA
      }
      m<-bigN-smalln
      Ea1<-(smalln*m/(bigN-1))
      Va1<-2*smalln*(smalln-1)*m*(m-1)/((bigN-3)*(bigN-1)*(bigN-1))
      dev<-(a1-Ea1)/sqrt(Va1)
      approx<-stats::pnorm(dev)
      list(a1=a1,Ea1=Ea1,Va1=Va1,dev=dev,pval=pval,approxpval=approx)
    }else{
      warning('Consider choosing a different covariate balance measure.')
    }
  }else if (!is.null(fdata)){
    stopifnot(!is.null(fz))
    stopifnot(all(fz%in%c(0,1)))

    if (is.vector(fdata)) fdata<-matrix(fdata,ncol=1)
    fd0<-fdata[which(fz==0),]
    fd1<-fdata[which(fz==1),]
    if (is.vector(fd0)) fd0<-matrix(fd0,length(fd0),1)
    if (is.vector(fd1)) fd1<-matrix(fd1,length(fd1),1)

    if (measure=='smd'){
      varf0<-apply(fd0,2,var,na.rm=TRUE)
      varf1<-apply(fd1,2,var,na.rm=TRUE)
      meanf0<-apply(fd0,2,mean,na.rm=TRUE)
      meanf1<-apply(fd1,2,mean,na.rm=TRUE)
      smdf<-(meanf1-meanf0)/sqrt((varf0+varf1)/2)

      r<-cbind(meanf1,meanf0,smdf)
      rownames(r)<-colnames(fdata)
      colnames(r)<-c('Treated Mean','Control Mean','SMD')
      r
    }else if (measure=='t.test'){
      nvar=ncol(fdata)
      meanf0<-apply(fd0,2,mean,na.rm=TRUE)
      meanf1<-apply(fd1,2,mean,na.rm=TRUE)
      t.test.pval.before=rep(0,nvar)
      for (i in 1:nvar){
        t.test.pval.before[i]=
          stats::t.test(fd1[,i],fd0[,i])$p.value
      }
      r=cbind(meanf1,meanf0,t.test.pval.before)
      rownames(r)<-colnames(fdata)
      colnames(r)<-c('Treated Mean','Control Mean','P-value')
      r
    }else if (measure=='crossmatch'){
      stop('The crossmatch option is only available for the matched data.')
    }else{
      warning('Consider choosing a different covariate balance measure.')
    }
  }else{
    stopifnot(!is.null(mz))
    stopifnot(all(mz%in%c(0,1)))

    if (is.vector(mdata)) mdata<-matrix(mdata,ncol=1)
    md0<-mdata[which(mz==0),]
    md1<-mdata[which(mz==1),]
    if (is.vector(md0)) md0<-matrix(md0,length(md0),1)
    if (is.vector(md1)) md1<-matrix(md1,length(md1),1)

    if (measure=='smd'){
      stop('fdata and fz are needed to calculate the SMDs.')
    }else if (measure=='t.test'){
      nvar=ncol(fdata)
      meanm0<-apply(md0,2,mean,na.rm=TRUE)
      meanm1<-apply(md1,2,mean,na.rm=TRUE)
      t.test.pval.after=rep(0,nvar)
      for (i in 1:nvar){
        t.test.pval.after[i]=stats::t.test(md1[,i],md0[,i])$p.value
      }
      r=cbind(meanm1,meanm0,t.test.pval.after)
      rownames(r)<-colnames(fdata)
      colnames(r)<-c('Treated Mean','Control Mean','P-value')
      r
    }else if (measure=='crossmatch'){
      D=maha.all(mdata)
      plainmatrix<-as.matrix(s.dist*max(1/min(D[D>0]),1)*D)
      diag(plainmatrix) <- 0
      mdm<-nbpMatching::distancematrix(plainmatrix)
      res<-nbpMatching::nonbimatch(mdm)

      mt_min<-pmin(as.numeric(res$matches$Group1.Row),
                   as.numeric(res$matches$Group2.Row))
      mt_max<-pmax(as.numeric(res$matches$Group1.Row),
                   as.numeric(res$matches$Group2.Row))
      z0<-mz[(mt_min>0)&(mt_max<=length(mz))]
      mt0<-factor(mt_min[(mt_min>0)&(mt_max<=length(mz))])
      tab<-table(factor(z0),mt0)
      a1<-sum(tab[1,]==1)

      crossmatchdist<-function(N,n){
        if (bigN%%2 == 1){
          stop("The number of subjects, bigN, should be even")
          return(NA)
        }
        I<-N/2
        cmd<-c()
        for (a1 in 0:I){
          a2 <- (n-a1)/2
          if ((floor(a2)==a2)&(a2>=0)){
            a0<-I-(a1+a2)
            if (a0>=0){
              pr<-factorial(I)/choose(N,n)
              pr<-pr*(2^a1)/(factorial(a0)*factorial(a1)*factorial(a2))
              cmd<-cbind(cmd,c(a0,a1,a2,pr))
            }
          }
        }
        cmd <- rbind(cmd,cumsum(cmd[4,]))
        cmd
      }

      bigN<-length(z0)
      smalln<-sum(z0)
      if (bigN<340){
        dist<-crossmatchdist(bigN,smalln)
        pval<-dist[5,dist[2,]==a1]
      }else{
        pval<-NA
      }
      m<-bigN-smalln
      Ea1<-(smalln*m/(bigN-1))
      Va1<-2*smalln*(smalln-1)*m*(m-1)/((bigN-3)*(bigN-1)*(bigN-1))
      dev<-(a1-Ea1)/sqrt(Va1)
      approx<-stats::pnorm(dev)
      list(a1=a1,Ea1=Ea1,Va1=Va1,dev=dev,pval=pval,approxpval=approx)
    }else{
      warning('Consider choosing a different covariate balance measure.')
    }
  }

}
