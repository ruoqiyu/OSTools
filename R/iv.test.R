iv.test <- function(y,treatment=NULL,instrument=NULL,mset=NULL,beta=0,
                    method=c("signedrank","ranksum"),Gamma=1,
                    alternative=c("two.sided","greater","less"),
                    exact=NULL,correct=TRUE,zero.method=c("wilcox","pratt"),
                    d=NULL,z=NULL) {
  method <- match.arg(method); alternative <- match.arg(alternative)
  zero.method <- match.arg(zero.method)
  .iv_scalar(beta,"beta")
  .iv_options(Gamma,exact,correct,zero.method,method)
  treatment.expr <- if (!is.null(d)) substitute(d) else substitute(treatment)
  inputs <- .iv_resolve_inputs(treatment,instrument,d,z)
  dat <- .iv_prepare(y,inputs$treatment,inputs$instrument,mset,method)
  ans <- .iv_evaluate(dat,beta,method,Gamma,alternative,exact,correct,zero.method)
  label <- if(method=="signedrank") "Wilcoxon signed-rank" else "Wilcoxon rank-sum"
  label <- paste(label,"test of IV-adjusted outcomes")
  if (Gamma>1) label <- paste(label,"with sensitivity bounds")
  else label <- paste(label,if(ans$exact) "(exact)" else "(normal approximation)")
  structure(list(statistic=ans$statistic,p.value=ans$p.value,
                 p.value.bounds=ans$bounds,null.value=c(beta=beta),
                 alternative=alternative,method=label,
                 data.name=paste(paste(deparse(substitute(y)),collapse=""),"adjusted by",
                                 paste(deparse(treatment.expr),collapse="")),
                 Gamma=Gamma,first.stage=dat$first.stage,n=dat$n,n.discarded=dat$excluded,
                 exact=ans$exact,correct=correct,zero.method=zero.method,call=match.call()),
            class=c("os_iv_test","htest"))
}

print.os_iv_test <- function(x,...) {
  NextMethod()
  cat("First-stage difference:",format(x$first.stage),"\n")
  if(x$Gamma>1) cat("Gamma:",format(x$Gamma),"; p-value bounds:",
                    paste(format(x$p.value.bounds),collapse=" to "),"\n")
  invisible(x)
}
