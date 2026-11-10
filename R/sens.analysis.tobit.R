sens.analysis.tobit <- function(y,z,mset,Gamma=1,tau=0,
    alternative=c("two-sided","greater","less"),method=c("mean","m"),
    inner=0,trim=3,lambda=1/2) {
  alternative <- match.arg(alternative);method <- match.arg(method)
  input <- .tobit.inputs(y,z,mset,Gamma,method,inner,trim,lambda)
  if (!is.numeric(tau) || length(tau)!=1L || !is.null(dim(tau)) ||
      !is.finite(tau) || tau<0 || tau>1) stop("tau must be a finite scalar in [0,1].")
  ans <- .tobit.evaluate(input,tau,alternative)
  ans$call <- match.call();class(ans) <- "os_tobit_test";ans
}

print.os_tobit_test <- function(x,...) {
  cat("Tobit-effect sensitivity test (positive tau is a reduction)\n")
  cat("method:",x$method," Gamma:",x$Gamma," tau:",x$tau,
      " alternative:",x$alternative,"\n")
  cat("p-value:",format.pval(x$p.value),"\n",x$approximation,"\n")
  invisible(x)
}
