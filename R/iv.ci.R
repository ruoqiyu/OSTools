.iv_grid_components <- function(grid, accepted) {
  if (!any(accepted)) return(matrix(numeric(),nrow=0L,ncol=2L,
                                   dimnames=list(NULL,c("lower","upper"))))
  starts <- which(accepted & !c(FALSE,utils::head(accepted,-1L)))
  ends <- which(accepted & !c(utils::tail(accepted,-1L),FALSE))
  cbind(lower=grid[starts],upper=grid[ends])
}

iv.ci <- function(y,treatment=NULL,instrument=NULL,mset=NULL,method=c("signedrank","ranksum"),
                  Gamma=1,level=0.95,grid=NULL,interval=c(-10,10),tol=0.01,
                  exact=NULL,correct=TRUE,zero.method=c("wilcox","pratt"),
                  d=NULL,z=NULL) {
  method <- match.arg(method); zero.method <- match.arg(zero.method)
  .iv_options(Gamma,exact,correct,zero.method,method)
  .iv_scalar(level,"level",lower=0,upper=1,lower.open=TRUE,upper.open=TRUE)
  if (is.null(grid)) {
    if (!is.numeric(interval) || length(interval)!=2L || any(!is.finite(interval)) || interval[1L]>=interval[2L])
      stop("interval must contain two increasing finite endpoints.",call.=FALSE)
    .iv_scalar(tol,"tol",lower=0,lower.open=TRUE)
    grid <- seq(interval[1L],interval[2L],by=tol)
  }
  if (!is.numeric(grid) || !is.null(dim(grid)) || length(grid)<2L ||
      any(!is.finite(grid)) || any(diff(grid)<=0))
    stop("grid must contain at least two strictly increasing finite values.",call.=FALSE)
  inputs <- .iv_resolve_inputs(treatment,instrument,d,z)
  dat <- .iv_prepare(y,inputs$treatment,inputs$instrument,mset,method)
  evaluate <- function(g) lapply(grid,function(b)
    .iv_evaluate(dat,b,method,g,"two.sided",exact,correct,zero.method,warn=FALSE))
  tests <- evaluate(Gamma)
  p <- vapply(tests,function(a) a$p.value,numeric(1))
  if (any(vapply(tests,function(a) a$fallback,logical(1))))
    warning("Some grid tests require a normal approximation because of ties, zero differences or Pratt scores.",call.=FALSE)
  accepted <- p>=(1-level)*(1-8*.Machine$double.eps)
  components <- .iv_grid_components(grid,accepted)
  boundary <- c(lower=accepted[1L],upper=utils::tail(accepted,1L))
  if (!any(accepted)) warning("No effects were accepted on the supplied grid; this does not establish an empty confidence set outside it.",call.=FALSE)
  if (any(boundary)) warning("The accepted set reaches a search boundary. Extend the grid; these grid limits are not confirmed finite confidence endpoints.",call.=FALSE)
  # Sensitivity sets retain the Gamma=1 grid point estimate, rather than
  # treating a plateau of worst-case p-values as a new causal point estimate.
  base.p <- if(Gamma==1) p else vapply(evaluate(1),function(a) a$p.value,numeric(1))
  estimate <- if(length(unique(base.p))==1L) NA_real_ else grid[which.max(base.p)]
  span <- if(nrow(components)) c(components[1L,1L],components[nrow(components),2L]) else c(NA_real_,NA_real_)
  structure(list(estimate=estimate,components=components,span=span,
                 touches.boundary=any(boundary),boundary.contact=boundary,
                 empty=!any(accepted),disconnected=nrow(components)>1L,
                 grid=grid,p.values=p,accepted=accepted,level=level,Gamma=Gamma,
                 first.stage=dat$first.stage,n=dat$n,n.discarded=dat$excluded,
                 method=method,exact=exact,correct=correct,zero.method=zero.method,
                 call=match.call()),class="os_iv_ci")
}

print.os_iv_ci <- function(x,...) {
  cat("Grid-based",format(100*x$level),"% IV confidence set\n",sep=" ")
  cat("Method:",x$method,"; Gamma:",format(x$Gamma),"\n")
  cat("Grid point estimate (Gamma = 1):",format(x$estimate),"\n")
  if(x$empty) cat("No accepted grid values.\n") else print(x$components,...)
  if(x$disconnected) cat("The accepted components are separate; their span fills rejected gaps.\n")
  cat("Touches search boundary:",x$touches.boundary,"\n")
  if(x$touches.boundary) cat("Extend the grid before interpreting its outer endpoints.\n")
  invisible(x)
}
