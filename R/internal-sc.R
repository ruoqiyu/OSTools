.sc_matrix <- function(x, name, min.units=1L) {
  if (is.data.frame(x)) {
    if (!all(vapply(x,is.numeric,logical(1L))))
      stop(name," must contain only numeric columns.",call.=FALSE)
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || nrow(x)<2L ||
      ncol(x)<min.units || any(!is.finite(x)))
    stop(name," must be a finite numeric matrix with at least two rows and ",
         min.units," column(s).",call.=FALSE)
  labels <- colnames(x)
  if (is.null(labels)) labels <- paste0("unit",seq_len(ncol(x)))
  if (anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels))
    stop(name," column names must be nonempty and unique.",call.=FALSE)
  colnames(x) <- labels
  x
}

.sc_index <- function(x, n, name) {
  if (is.logical(x) && is.null(dim(x)) && length(x)==n && !anyNA(x)) {
    ans <- which(x)
  } else if (is.numeric(x) && is.null(dim(x)) && length(x)>0L &&
             all(is.finite(x)) && all(x==floor(x)) &&
             all(x>=1 & x<=n) && !anyDuplicated(x)) {
    ans <- sort(as.integer(x))
  } else stop(name," must be a logical vector of length nrow(donors/outcomes) or unique row indices.",call.=FALSE)
  if (!length(ans)) stop(name," must select at least one period.",call.=FALSE)
  ans
}

.sc_periods <- function(pre, post, n) {
  pre <- .sc_index(pre,n,"pre")
  if (is.null(post)) post <- setdiff(seq_len(n),pre)
  post <- .sc_index(post,n,"post")
  if (length(intersect(pre,post)))
    stop("pre and post must not overlap.",call.=FALSE)
  if (max(pre)>=min(post))
    stop("Rows must be in chronological order, with every pre period before every post period.",call.=FALSE)
  list(pre=pre,post=post)
}

.sc_ridge <- function(ridge) {
  if (!is.numeric(ridge) || !is.null(dim(ridge)) || length(ridge)!=1L ||
      !is.finite(ridge) || ridge<0)
    stop("ridge must be one finite nonnegative number.",call.=FALSE)
}

.sc_weights <- function(y, donors, ridge) {
  j <- ncol(donors)
  D <- crossprod(donors)+diag(ridge,j)
  d <- crossprod(donors,y)
  if (any(!is.finite(D)) || any(!is.finite(d)))
    stop("The quadratic program overflowed; rescale the outcomes.",call.=FALSE)
  ans <- tryCatch(quadprog::solve.QP(D,d,cbind(rep(1,j),diag(j)),
                                    c(1,rep(0,j)),meq=1),
                  error=function(e) stop("Synthetic-control optimization failed: ",
                    conditionMessage(e)," Rescale the outcomes or use a larger ridge.",call.=FALSE))
  w <- ans$solution
  if (any(!is.finite(w)) || min(w)< -1e-7 || abs(sum(w)-1)>1e-7)
    stop("The optimizer did not return feasible donor weights.",call.=FALSE)
  # Remove tiny numerical violations, using the chapter's normalization.
  w <- pmax(w,0)
  stats::setNames(w/sum(w),colnames(donors))
}
