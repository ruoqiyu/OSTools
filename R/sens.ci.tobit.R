sens.ci.tobit <- function(y,z,mset,Gamma=1,alpha=0.05,
    alternative=c("two-sided","greater","less"),method=c("mean","m"),
    inner=0,trim=3,lambda=1/2,grid=seq(0,1,length.out=201),tol=1e-6) {
  alternative <- match.arg(alternative);method <- match.arg(method)
  input <- .tobit.inputs(y,z,mset,Gamma,method,inner,trim,lambda)
  scalar <- function(v) is.numeric(v) && length(v)==1L && is.null(dim(v)) && is.finite(v)
  if (!scalar(alpha) || alpha<=0 || alpha>=1) stop("alpha must lie strictly between zero and one.")
  if (!scalar(tol) || tol<=0 || tol>0.01) stop("tol must lie in (0,0.01].")
  if (!is.numeric(grid) || !is.null(dim(grid)) || length(grid)<2L ||
      any(!is.finite(grid)) || any(grid<0 | grid>1)) stop("grid must contain at least two finite values in [0,1].")
  grid <- sort(unique(c(0,grid,1)))
  cache <- new.env(parent=emptyenv())
  evaluate <- function(t) {
    key <- sprintf("%.17g",t)
    if (!exists(key,cache,inherits=FALSE)) {
      a <- .tobit.evaluate(input,t,alternative)
      assign(key,c(a$p.value-alpha,-a$scores),cache)
    }
    get(key,cache,inherits=FALSE)
  }
  confidence <- .tobit.regions(evaluate,grid,1L,tol,1e-12)
  # Point-estimate ranges solve the sensitivity estimating inequalities. At
  # Gamma=1 these reduce to a zero of the centered statistic; flat regions stay.
  score.scale <- max(1,max(abs(vapply(grid,function(t) evaluate(t)[2:3],numeric(2)))))
  estimates <- .tobit.regions(evaluate,grid,2:3,tol,1e-8*score.scale)
  envelope <- function(x) if (!nrow(x)) c(lower=NA_real_,upper=NA_real_) else
    c(lower=min(x[,1]),upper=max(x[,2]))
  table <- data.frame(tau=grid,p.value=vapply(grid,function(t) evaluate(t)[1]+alpha,numeric(1)))
  ans <- list(confidence.set=confidence,conf.int=envelope(confidence),
    estimate.set=estimates,point.estimates=envelope(estimates),
    Gamma=Gamma,alpha=alpha,alternative=alternative,method=method,
    grid=table,grid.resolution=max(diff(grid)),tol=tol,numerical=TRUE,
    approximation="Numerical test inversion using the separable Normal approximation; scan resolution matters.",
    call=match.call())
  class(ans) <- "os_tobit_ci";ans
}

print.os_tobit_ci <- function(x,...) {
  cat("Numerical Tobit-effect confidence set\n")
  cat("Gamma:",x$Gamma," confidence level:",1-x$alpha," method:",x$method,"\n")
  print(x$confidence.set)
  cat("Point-estimate set:\n");print(x$estimate.set)
  cat("Largest scan spacing:",x$grid.resolution,"; root tolerance:",x$tol,"\n")
  invisible(x)
}
