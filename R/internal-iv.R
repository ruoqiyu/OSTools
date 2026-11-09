.iv_resolve_inputs <- function(treatment, instrument, d, z) {
  # The legacy argument z remains the instrument. Never reinterpret it as the
  # treatment: explicit role names avoid an ambiguous D/Z -> Z/A migration.
  if (!is.null(treatment) && !is.null(d))
    stop("Supply treatment or the legacy d argument, not both.",call.=FALSE)
  if (!is.null(instrument) && !is.null(z))
    stop("Supply instrument or the legacy z argument, not both.",call.=FALSE)
  if (is.null(treatment)) treatment <- d
  if (is.null(instrument)) instrument <- z
  if (is.null(treatment))
    stop("Supply treatment (Z); d is retained as a legacy alias.",call.=FALSE)
  list(treatment=treatment,instrument=instrument)
}

.iv_scalar <- function(x, name, lower=-Inf, upper=Inf, lower.open=FALSE, upper.open=FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x)!=1L || !is.finite(x) ||
      (if (lower.open) x<=lower else x<lower) ||
      (if (upper.open) x>=upper else x>upper))
    stop(name, " must be one finite number in the permitted range.",call.=FALSE)
  invisible(x)
}

.iv_logical <- function(x, name, nullable=FALSE) {
  if (nullable && is.null(x)) return(invisible(x))
  if (!is.logical(x) || length(x)!=1L || is.na(x))
    stop(name, " must be TRUE or FALSE",if(nullable) " (or NULL)" else "", ".",call.=FALSE)
  invisible(x)
}

.iv_vector <- function(x, name) {
  if (!(is.numeric(x) || is.logical(x)) || !is.null(dim(x)) || !length(x) || any(!is.finite(x)))
    stop(name," must be a nonempty finite numeric vector.",call.=FALSE)
}

.iv_prepare <- function(y, treatment, instrument, mset, method) {
  excluded <- 0L
  if (is.logical(y)) storage.mode(y) <- "double"
  if (is.logical(treatment)) storage.mode(treatment) <- "double"
  if (is.matrix(y) || is.matrix(treatment)) {
    if (method!="signedrank" || !is.null(instrument) || !is.null(mset) ||
        !is.matrix(y) || !is.matrix(treatment) || !is.numeric(y) || !is.numeric(treatment) ||
        !identical(dim(y),dim(treatment)) || ncol(y)!=2L || nrow(y)<1L ||
        any(!is.finite(y)) || any(!is.finite(treatment)))
      stop("Paired matrix inputs must be finite numeric y and treatment matrices with identical dimensions and two columns; omit instrument and mset.",call.=FALSE)
    y1 <- y[,1L]; y0 <- y[,2L]; t1 <- treatment[,1L]; t0 <- treatment[,2L]
  } else {
    .iv_vector(y,"y"); .iv_vector(treatment,"treatment")
    n <- length(y)
    if (length(treatment)!=n || !(is.numeric(instrument) || is.logical(instrument)) || !is.null(dim(instrument)) ||
        length(instrument)!=n || any(!is.finite(instrument)))
      stop("Vector inputs need equally long y, treatment and finite numeric or logical instrument.",call.=FALSE)
    if (method=="ranksum") {
      if (!is.null(mset)) stop("mset is only supported for signedrank paired analysis.",call.=FALSE)
      levels <- sort(unique(instrument))
      if (length(levels)!=2L) stop("A ranksum analysis requires exactly two instrument levels.",call.=FALSE)
      hi <- which(instrument==levels[2L]); lo <- which(instrument==levels[1L])
    } else {
      if (!(is.atomic(mset) || is.factor(mset)) || !is.null(dim(mset)) || length(mset)!=n)
        stop("A signedrank vector analysis requires mset labels for matched pairs.",call.=FALSE)
      keep <- !is.na(mset); excluded <- sum(!keep)
      if (!any(keep)) stop("mset contains no matched subjects.",call.=FALSE)
      ids <- as.character(mset[keep]); rows <- which(keep)
      sets <- lapply(unique(ids),function(id) rows[ids==id])
      if (any(lengths(sets)!=2L)) stop("Every retained mset must contain exactly two subjects.",call.=FALSE)
      if (any(vapply(sets,function(i) instrument[i[1L]]==instrument[i[2L]],logical(1))))
        stop("The instrument must differ within every retained pair.",call.=FALSE)
      hi <- vapply(sets,function(i) i[which.max(instrument[i])],integer(1))
      lo <- vapply(sets,function(i) i[which.min(instrument[i])],integer(1))
    }
    y1 <- y[hi]; y0 <- y[lo]; t1 <- treatment[hi]; t0 <- treatment[lo]
  }
  list(y1=y1,y0=y0,t1=t1,t0=t0,excluded=excluded,
       first.stage=mean(t1)-mean(t0),n=c(higher=length(y1),lower=length(y0)))
}

.iv_options <- function(Gamma, exact, correct, zero.method, method) {
  .iv_scalar(Gamma,"Gamma",lower=1)
  .iv_logical(exact,"exact",nullable=TRUE); .iv_logical(correct,"correct")
  if (Gamma>1 && method!="signedrank")
    stop("Gamma greater than 1 is supported only for paired signedrank analysis.",call.=FALSE)
  if (Gamma>1 && isTRUE(exact))
    stop("Sensitivity bounds use a normal approximation; exact=TRUE requires Gamma=1.",call.=FALSE)
  invisible(NULL)
}

.iv_rank_bounds <- function(difference, Gamma, alternative, correct, zero.method) {
  if (zero.method=="wilcox") difference <- difference[difference!=0]
  if (!length(difference) || all(difference==0))
    return(list(p.value=1,bounds=c(lower=1,upper=1),statistic=c(V=0),exact=FALSE))
  scores <- rank(abs(difference)); active <- difference!=0
  W <- sum(scores[difference>0]); score.sum <- sum(scores[active])
  q <- 1/(1+Gamma)
  low <- q*score.sum; high <- (1-q)*score.sum
  se <- sqrt(q*(1-q)*sum(scores[active]^2))
  cc <- if (correct) 0.5 else 0
  if (alternative=="greater") {
    bounds <- stats::pnorm((W-c(low,high)-cc)/se,lower.tail=FALSE)
  } else if (alternative=="less") {
    bounds <- stats::pnorm((W-c(high,low)+cc)/se)
  } else {
    far <- max(abs(W-c(low,high)))
    near <- if(W>=low && W<=high) 0 else min(abs(W-c(low,high)))
    bounds <- 2*stats::pnorm(-pmax(c(far,near)-cc,0)/se)
  }
  bounds <- stats::setNames(pmin(1,pmax(0,bounds)),c("lower","upper"))
  list(p.value=unname(bounds[2L]),bounds=bounds,statistic=c(V=W),exact=FALSE)
}

.iv_evaluate <- function(dat, beta, method, Gamma, alternative, exact, correct, zero.method,
                         warn=TRUE) {
  r1 <- dat$y1-beta*dat$t1; r0 <- dat$y0-beta*dat$t0
  if (any(!is.finite(c(r1,r0))))
    stop("The adjusted outcomes overflowed; rescale the data or effect grid.",call.=FALSE)
  if (method=="signedrank") {
    difference <- r1-r0
    if (any(!is.finite(difference))) stop("Paired differences overflowed; rescale the data.",call.=FALSE)
    if (all(difference==0))
      return(list(p.value=1,bounds=c(lower=1,upper=1),statistic=c(V=0),exact=FALSE,fallback=FALSE))
    untied <- !anyDuplicated(abs(difference)) && !any(difference==0)
    eligible <- zero.method=="wilcox" && untied
    auto <- sum(difference!=0)<50L && eligible
  } else {
    eligible <- !anyDuplicated(c(r1,r0))
    auto <- length(r1)<50L && length(r0)<50L && eligible
  }
  fallback <- isTRUE(exact) && !eligible
  if (fallback && warn)
    warning("Exact rank inference is not implemented here for ties, zero differences or Pratt scores; using a normal approximation.",call.=FALSE)
  exact.used <- Gamma==1 && eligible && if(is.null(exact)) auto else isTRUE(exact)
  if (method=="signedrank" && (Gamma>1 || zero.method=="pratt")) {
    ans <- .iv_rank_bounds(difference,Gamma,alternative,correct,zero.method)
  } else {
    test <- stats::wilcox.test(r1,r0,paired=method=="signedrank",alternative=alternative,
                              exact=exact.used,correct=correct)
    p <- test$p.value
    # Identical tied outcomes have no evidence against the null, even when the
    # zero-variance normal approximation in wilcox.test returns NaN.
    if (is.nan(p) && length(unique(c(r1,r0)))==1L) p <- 1
    if (!is.finite(p)) stop("The rank test did not produce a finite p-value.",call.=FALSE)
    ans <- list(p.value=p,bounds=c(lower=p,upper=p),statistic=test$statistic,exact=exact.used)
  }
  ans$fallback <- fallback
  ans
}
