.os_pvalues <- function(p, log.p) {
  if (!is.logical(log.p) || length(log.p) != 1L || is.na(log.p))
    stop("log.p must be TRUE or FALSE.")
  if (!is.numeric(p) || !is.null(dim(p)) || !length(p) || anyNA(p))
    stop("p must be a nonempty numeric vector without missing values.")
  if (log.p) {
    if (any(p > 0)) stop("Log p-values must be nonpositive (negative infinity is allowed).")
    lp <- p
    probabilities <- exp(p)
  } else {
    if (any(!is.finite(p)) || any(p < 0 | p > 1))
      stop("p-values must be finite and between zero and one.")
    lp <- log(p)
    probabilities <- p
  }
  labels <- names(p)
  if (is.null(labels)) labels <- paste0("H", seq_along(p))
  if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))
    stop("p-value names must be unique and nonempty.")
  names(probabilities) <- names(lp) <- labels
  list(p = probabilities, log.p = lp, labels = labels)
}

.os_combination_method <- function(method, trunc) {
  method <- match.arg(method, c("fisher", "truncated"))
  if (!is.numeric(trunc) || length(trunc) != 1L || !is.finite(trunc) ||
      trunc <= 0 || trunc > 1)
    stop("trunc must be a finite number strictly above zero and at most one.")
  list(method = method, cutoff = if (method == "fisher") 1 else trunc)
}

.os_combine_logp <- function(lp, cutoff) {
  n <- length(lp)
  log.cutoff <- log(cutoff)
  log.product <- sum(lp[lp <= log.cutoff])
  if (log.product == 0) {
    # Includes the empty product and the point mass at product = 1.
    log.p <- 0
  } else if (n == 1L) {
    log.p <- log.product
  } else if (is.infinite(log.product)) {
    log.p <- -Inf
  } else if (cutoff == 1) {
    # Gamma(n, 1) is half a chi-squared variable with 2*n degrees of freedom.
    # Its argument remains finite even when -2*log.product would overflow.
    log.p <- stats::pgamma(-log.product, shape = n,
                          lower.tail = FALSE, log.p = TRUE)
  } else {
    # Conditional on J retained uniforms, W = cutoff^J * product(V_1,...,V_J).
    # Sum nonnegative binomial-mixture terms on the log scale.
    j <- seq_len(n)
    threshold <- pmax(0, j * log.cutoff - log.product)
    terms <- stats::dbinom(j, size = n, prob = cutoff, log = TRUE) +
      stats::pgamma(threshold, shape = j, lower.tail = FALSE, log.p = TRUE)
    largest <- max(terms)
    log.p <- if (is.infinite(largest)) -Inf else
      min(0, largest + log(sum(exp(terms - largest))))
  }
  list(log.p.value = log.p, log.product = log.product,
       retained = sum(lp <= log.cutoff))
}
