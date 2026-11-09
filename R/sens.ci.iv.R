sens.ci.iv <- function(encouraged.treatment,unencouraged.treatment,
                       encouraged.response,unencouraged.response,
                       Gamma=1,alpha=0.05,tol=0.01,interval=c(-10,10),details=FALSE) {
  inputs <- list(encouraged.treatment,unencouraged.treatment,
                 encouraged.response,unencouraged.response)
  for(i in seq_along(inputs)) .iv_vector(inputs[[i]],names(formals(sens.ci.iv))[i])
  if(length(unique(lengths(inputs)))!=1L)
    stop("All four paired input vectors must have the same length.",call.=FALSE)
  .iv_scalar(alpha,"alpha",lower=0,upper=1,lower.open=TRUE,upper.open=TRUE)
  .iv_logical(details,"details")
  ans <- iv.ci(y=cbind(encouraged.response,unencouraged.response),
               treatment=cbind(encouraged.treatment,unencouraged.treatment),
               Gamma=Gamma,level=1-alpha,interval=interval,tol=tol,
               method="signedrank",exact=FALSE,correct=FALSE,zero.method="pratt")
  if(details) return(ans)
  if(ans$disconnected)
    warning("The accepted set is disconnected; the legacy endpoint vector gives its span. Use details=TRUE for all components.",call.=FALSE)
  ans$span
}
