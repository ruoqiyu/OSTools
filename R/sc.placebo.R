sc.placebo <- function(outcomes, treated, pre, post=NULL, ridge=1e-8) {
  call <- match.call()
  outcomes <- .sc_matrix(outcomes,"outcomes",min.units=2L)
  labels <- colnames(outcomes)
  if (is.character(treated) && length(treated)==1L && !is.na(treated)) {
    treated <- match(treated,labels)
  }
  if (!is.numeric(treated) || !is.null(dim(treated)) || length(treated)!=1L ||
      !is.finite(treated) || treated!=floor(treated) ||
      treated<1 || treated>ncol(outcomes))
    stop("treated must identify one outcome column by its name or index.",call.=FALSE)
  treated <- as.integer(treated)
  periods <- .sc_periods(pre,post,nrow(outcomes))
  .sc_ridge(ridge)
  n <- ncol(outcomes)
  weights <- matrix(0,n,n,dimnames=list(labels,labels))
  gaps <- matrix(NA_real_,nrow(outcomes),n,dimnames=dimnames(outcomes))
  fits <- lapply(seq_len(n),function(k) {
    fit <- sc.fit(outcomes[,k],outcomes[,-k,drop=FALSE],
                  periods$pre,periods$post,ridge=ridge)
    weights[-k,k] <<- fit$weights
    gaps[,k] <<- fit$gap
    fit
  })
  before <- vapply(fits,function(x) x$pre.rmspe,numeric(1L))
  after <- vapply(fits,function(x) x$post.rmspe,numeric(1L))
  ratios <- after/before
  # A zero / zero ratio is undefined; never silently drop that candidate.
  ratios[before==0 & after==0] <- NA_real_
  names(ratios) <- labels
  defined <- !anyNA(ratios)
  rank <- if (defined) as.integer(sum(ratios>=ratios[treated])) else NA_integer_
  tab <- data.frame(unit=labels,pre.rmspe=before,post.rmspe=after,
                    ratio=unname(ratios),
                    post.mean.gap=vapply(fits,function(x) x$post.mean.gap,numeric(1L)),
                    stringsAsFactors=FALSE)
  structure(list(treated=labels[treated],statistics=tab,ratios=ratios,
                 rank=rank,rank.fraction=rank/n,reference.size=n,
                 rank.defined=defined,weights=weights,gaps=gaps,
                 pre=periods$pre,post=periods$post,ridge=ridge,
                 interpretation="Descriptive rank fraction; not automatically a randomization p-value.",
                 call=call),class="os_sc_placebo")
}

print.os_sc_placebo <- function(x,...) {
  cat("Descriptive in-space synthetic-control placebos\n")
  cat("Treated unit:",x$treated,"; reference units:",x$reference.size,"\n")
  print(c(rank=x$rank,rank.fraction=x$rank.fraction),...)
  if (!x$rank.defined)
    cat("Rank undefined: a candidate has zero pre and post RMSPE; no units were dropped.\n")
  cat(x$interpretation,"\n")
  invisible(x)
}
