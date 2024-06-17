add.caliper<-function(dist,z,dx,rg, stdev = FALSE, penalty = 1000, constant = TRUE){

  if (length(rg)==1) rg=c(-abs(rg),abs(rg))
  stopifnot(is.vector(rg)&(length(rg)==2))
  stopifnot((rg[1]<=0)&(rg[2]>=0))

  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))

  dx1=dx[z==1]
  dx0=dx[z==0]
  # Standardize p using an equally weighted pooled variance
  v1=stats::var(dx1)
  v2=stats::var(dx0)
  sp=sqrt((v1+v2)/2)
  stopifnot(sp>0)

  if (stdev) rg=rg *sp

  m=sum(z)

  if (is.matrix(dist)){
    dif=outer(dx1,dx0,"-")
    if (constant) dif=(dif>rg[2])+(dif<rg[1])
    else dif=(dif-rg[2])*(dif>rg[2])+(rg[1]-dif)*(dif<rg[1])
    dist=dist+dif*penalty
    return (dist)
  }else if (is.list(dist)){
    dif=dx1[dist$start]-dx0[dist$end-m]
    if (constant) dif=(dif>rg[2])+(dif>rg[2])
    else dif=(dif-rg[2])*(dif>rg[2])+(rg[1]-dif)*(dif<rg[1])
    d=dist$d+dif*penalty
    return (list(d=d,start=dist$start,end=dist$end))
  }else{
    warning("dist has to be a matrix or a list.")
  }
}
