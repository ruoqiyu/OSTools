match.os<-function(z,dist,dat,p=rep(1,length(z)),exact=NULL,fine=rep(1,length(z)),ncontrol=1,
                   penalty=ifelse(is.matrix(dist),round(max(dist)*1000),round(max(dist$d)*1000)),
                   s.cost=100,subX=NULL,solver='rlemon'){
  #Check input
  stopifnot(is.data.frame(dat))
  stopifnot(is.vector(z))
  if (!is.numeric(fine)){
    if (!is.factor(fine)) fine=as.factor(fine)
    levels(fine)<-1:nlevels(fine)
    fine<-as.numeric(fine)
  }

  stopifnot(all((z==1)|(z==0)))
  stopifnot((ncontrol==round(ncontrol))&(ncontrol>=1))
  n<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(ncontr>=(ncontrol*ntreat))
  stopifnot(length(z)==length(fine))

  stopifnot(length(z)==(dim(dat)[1]))

  if (is.matrix(dist)){
    distance<-t(dist)
    dim(distance)<-c(1,ntreat*ncontr)
    distance<-as.vector(distance)
    start<-rep(1:ntreat,each=ncontr)
    end<-rep((ntreat+1):n,ntreat)
    d0<-distance
    distance<-distance[which(d0<Inf)]
    start<-start[which(d0<Inf)]
    end<-end[which(d0<Inf)]
    dist=list(d=distance,start=start,end=end)
  }

  if (!is.null(subX)){
    if (is.factor(subX)){
      levels(subX)<-1:nlevels(subX)
      subX<-as.integer(subX)
    }
    stopifnot(is.vector(subX))
  }

  #sort input
  if (is.null(exact)){
    o<-order(1-p)
  }else{
    o<-order(exact,1-p)
    exact<-exact[o]
  }

  z<-z[o]
  p<-p[o]
  fine<-fine[o]
  dat<-dat[o,]

  #Must have treated first
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    o<-order(1-z)
    z<-z[o]
    dat<-dat[o,]
    fine<-fine[o]
    if (!is.null(subX)) subX<-subX[o]
  }

  #do match
  net<-function(z,dist,ncontrol=1,fine=rep(1,length(z)),
                penalty=ifelse(is.matrix(dist),round(max(dist)*1000),round(max(dist$d)*1000)),
                s.cost=100,subX=NULL){

    #check input
    stopifnot(is.vector(z))
    stopifnot(is.vector(fine))
    fine<-as.numeric(fine)
    stopifnot(all((z==1)|(z==0)))
    stopifnot((ncontrol==round(ncontrol))&(ncontrol>=1))
    nobs<-length(z)
    ntreat<-sum(z)
    ncontr<-sum(1-z)
    stopifnot(ncontr>=(ncontrol*ntreat))
    stopifnot(nobs==length(fine))

    if (!is.null(subX)){
      if (is.factor(subX)){
        levels(subX)<-1:nlevels(subX)
        subX<-as.integer(subX)
      }
      stopifnot(is.vector(subX))
    }

    #create basic treated-vs-control bipartite graph
    fine1=fine[z==1]
    fine0=fine[z==0]
    if (!is.null(subX)){
      subX1<-subX[z==1]
      subX0<-subX[z==0]
    }
    startn<-dist$start
    endn<-dist$end
    cost<-dist$d
    tcarcs<-length(startn) # number of treatment-control arcs
    ucap<-rep(1,tcarcs)
    b<-rep(ncontrol,ntreat) #supply for treated nodes
    b<-c(b,rep(0,ncontr)) #flow conservation at control nodes
    #Make costs integer
    if (any(cost<0)) cost<-cost-min(cost)
    #cost<-round(100*cost)
    cost<-round(cost*s.cost)

    #create a duplicate for each control to make sure each control is only used once
    startn=c(startn,(ntreat+1):nobs)
    endn=c(endn,(nobs+1):(nobs+ncontr))
    cost=c(cost,rep(0,ncontr))
    ucap=c(ucap,rep(1,ncontr))
    b<-c(b,rep(0,ncontr))

    #Add a node to take extras of subset for fine balance category k
    if (!is.null(subX)){
      tbs<-table(z,subX)
      ncs<-as.vector(tbs[1,])
      nts<-as.vector(tbs[2,])
      nwants<-nts*ncontrol #desired number
      nlows<-pmin(ncs,nwants) #available number
      nextras<-nwants-nlows #gap between desired and available
      sublevels<-as.vector(as.numeric(colnames(tbs)))
      extras<-c()
      for (kk in 1:length(nlows)){
        if (nextras[kk]>0){
          extrak<-length(b)+1
          extras<-c(extras,extrak)
          b<-c(b,0)
          who<-subX1==sublevels[kk]
          if (sum(who)>0){
            startn<-c(startn,which(who))
            endn<-c(endn,rep(extrak,sum(who)))
            ucap<-c(ucap,rep(1,sum(who)))
            cost<-c(cost,rep(0,sum(who)))
          }
        }
      }
    }

    #Add structure to the bipartite graph for near fine balance

    tb<-table(z,fine)
    nc<-as.vector(tb[1,])
    nt<-as.vector(tb[2,])
    nwant<-nt*ncontrol #desired number
    nlow<-pmin(nc,nwant) #available number
    nextra<-sum(nwant-nlow) #gap between desired and available
    finelevels<-as.vector(as.numeric(colnames(tb)))

    #Add a node for fine balance category k
    sinks<-NULL
    for (k in 1:length(nlow)){
      if (nlow[k]>0){
        sinkk<-length(b)+1
        sinks<-c(sinks,sinkk)
        who0<-fine0==finelevels[k]
        b<-c(b,0)
        if (sum(who0)>0){
          startn<-c(startn,rep(nobs,sum(who0))+which(who0))
          endn<-c(endn,rep(sinkk,sum(who0)))
          ucap<-c(ucap,rep(1,sum(who0)))
          cost<-c(cost,rep(0,sum(who0)))
        }
      }
    }

    #Add a node to take the extras
    sinkex<-length(b)+1
    b<-c(b,0)
    startn<-c(startn,(nobs+1):(nobs+ncontr))
    endn<-c(endn,rep(sinkex,ncontr))
    ucap<-c(ucap,rep(1,ncontr))
    cost<-c(cost,rep(0,ncontr))

    #Add a sink
    finalsink<-length(b)+1
    b<-c(b,-ntreat*ncontrol) #finalsink absorbs all flow
    #Connect balance nodes to finalsink
    startn<-c(startn,sinks)
    endn<-c(endn,rep(finalsink,length(sinks)))
    ucap<-c(ucap,nlow[nlow>0])
    cost<-c(cost,rep(0,length(sinks)))

    if(!is.null(subX)){
      startn<-c(startn,extras)
      endn<-c(endn,rep(finalsink,length(extras)))
      ucap<-c(ucap,nextras[nextras>0])
      cost<-c(cost,rep(0,length(extras)))
    }

    #Connect sinkex to finalsink
    startn<-c(startn,sinkex)
    endn<-c(endn,finalsink)
    ucap<-c(ucap,ntreat*ncontrol)
    #ucap<-c(ucap,nextra)
    cost<-c(cost,penalty)

    #Make costs integer
    cost<-round(cost)
    net<-list(startn=startn,endn=endn,ucap=ucap,b=b,cost=cost,tcarcs=tcarcs)
    net
  }

  net<-net(z,dist,ncontrol,fine,penalty,s.cost,subX)
  if (any(net$cost==Inf)) net$cost[net$cost==Inf]<-2*max(net$cost[net$cost!=Inf])

  callrelax <- function (net, solver = 'rlemon'){
    startn <- net$startn
    endn <- net$endn
    ucap <- net$ucap
    b <- net$b
    cost <- net$cost
    stopifnot(length(startn) == length(endn))
    stopifnot(length(startn) == length(ucap))
    stopifnot(length(startn) == length(cost))
    stopifnot(min(c(startn, endn)) >= 1)
    stopifnot(max(c(startn, endn)) <= length(b))
    stopifnot(all(startn != endn))

    nnodes <- length(b)
    if(solver == 'rrelaxiv'){
      if(requireNamespace('rrelaxiv', quietly = TRUE)) {
        rout <- rrelaxiv::RELAX_IV(startnodes = as.integer(startn),
                                   endnodes = as.integer(endn),
                                   arccosts = as.integer(cost),
                                   arccapacity = as.integer(ucap),
                                   supply = as.integer(b))
        return.obj <- list(crash = 0,
                           feasible = !all(rout == 0),
                           x = rout)
        return(return.obj)
      } else {
        solver = 'rlemon'
        warning('Package rrelaxiv not available, using rlemon instead.')
      }
    }
    if(solver == 'rlemon'){
      lout <- rlemon::MinCostFlow(arcSources = as.integer(startn),
                                  arcTargets = as.integer(endn),
                                  arcCapacities = as.integer(ucap),
                                  arcCosts = as.integer(cost),
                                  nodeSupplies = as.integer(b),
                                  numNodes = max(c(startn, endn)),
                                  algorithm = 'CycleCancelling')
      return.obj <- list(crash = 0,
                         feasible = !all(lout[[1]] == 0),
                         x = lout[[1]])
      return(return.obj)
    }else{
      stop(
        'Argument to solver not recognized: please use one of rlemon and rrelaxiv')
    }
  }

  output<-callrelax(net,solver)

  if (output$feasible!=1){
    warning("Match is infeasible.  Change dist or ncontrol to obtain a feasible match.")
    m<-list(feasible=output$feasible,d=NULL)
  }else{
    x<-output$x[1:net$tcarcs]
    treated<-net$startn[1:net$tcarcs]
    control<-net$endn[1:net$tcarcs]
    treated=treated[which(x==1)]
    control=control[which(x==1)]
    match.df=data.frame('treat'=treated,'control'=control)
    match.df$treat<-as.factor(as.character(match.df$treat))
    matches<-as.matrix(plyr::daply(match.df, plyr::.(match.df$treat),
                                   function(treat.edges) treat.edges$control,.drop_o=FALSE))
    id1<-(1:n)[z==1]
    id0<-(1:n)[z==0]
    matchid<-matrix(c(id1[as.numeric(row.names(matches))],id0[as.vector((matches-sum(z)))]),
                    ncol=ncontrol+1)
    matchid<-as.vector(t(matchid))
    dat1<-dat[matchid,]
    zm<-z[matchid]
    mset<-rep(1:nrow(matches),each=ncontrol+1)
    dat1$mset<-mset
    #dat1<-cbind(mset,dat1)
    m<-list(feasible=output$feasible,data=dat1,x=x)
  }

  if(m[[1]]==0) {
    warning("The match you requested is infeasible, reconsider caliper or ncontrol or exact for distance")
  }
  m
}
