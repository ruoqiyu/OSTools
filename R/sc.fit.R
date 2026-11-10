sc.fit <- function(y, donors, pre, post=NULL, ridge=1e-8) {
  call <- match.call()
  donors <- .sc_matrix(donors,"donors")
  if (!is.numeric(y) || !is.null(dim(y)) || length(y)!=nrow(donors) ||
      any(!is.finite(y)))
    stop("y must be a finite numeric vector with one value per donor row.",call.=FALSE)
  periods <- .sc_periods(pre,post,length(y))
  .sc_ridge(ridge)
  w <- .sc_weights(y[periods$pre],donors[periods$pre,,drop=FALSE],ridge)
  synthetic <- as.vector(donors %*% w)
  names(synthetic) <- names(y)
  gap <- y-synthetic
  if (any(!is.finite(gap)))
    stop("Predictions overflowed; rescale the outcomes.",call.=FALSE)
  pre.rmspe <- sqrt(mean(gap[periods$pre]^2))
  post.rmspe <- sqrt(mean(gap[periods$post]^2))
  if (!is.finite(pre.rmspe) || !is.finite(post.rmspe))
    stop("RMSPE calculation overflowed; rescale the outcomes.",call.=FALSE)
  structure(list(weights=w,synthetic=synthetic,gap=gap,outcome=y,
                 pre=periods$pre,post=periods$post,
                 pre.rmspe=pre.rmspe,post.rmspe=post.rmspe,
                 post.mean.gap=mean(gap[periods$post]),
                 ridge=ridge,call=call),class="os_sc_fit")
}

print.os_sc_fit <- function(x,...) {
  cat("Outcome-only synthetic control\n")
  cat("Donors:",length(x$weights),"; fitting periods:",length(x$pre),
      "; evaluation periods:",length(x$post),"\n")
  print(c(pre.rmspe=x$pre.rmspe,post.rmspe=x$post.rmspe,
          post.mean.gap=x$post.mean.gap),...)
  invisible(x)
}
