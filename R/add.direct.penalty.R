add.direct.penalty<-function(dist,z,dx,positive = TRUE,penalty=1){

  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))

  m <- sum(z)
  dx0 <- dx[z==0]
  dx1 <- dx[z==1]
  
  if (is.matrix(dist)){
    dif=outer(dx1,dx0,"-")
    if (positive) dist<-dist+(dif>0)*penalty
    else dist<-dist+(dif<0)*penalty
    return (dist)
  }else if (is.list(dist)){
    dif<-dx1[dist$start]-dx0[dist$end-m]
    if (positive) d<-dist$d+as.numeric(dif>0)*penalty
    else d<-dist$d+as.numeric(dif<0)*penalty
    return (list(d=d,start=dist$start,end=dist$end))
  }else{
    warning("dist has to be a matrix or a list.")
  }
}
