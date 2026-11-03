partial.conjunction <- function(p, method = c("fisher", "truncated"),
                                trunc = 0.2, alpha = 0.05,
                                r.followup = NULL, log.p = FALSE) {
  method <- match.arg(method)
  settings <- .os_combination_method(method, trunc)
  inputs <- .os_pvalues(p, log.p)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1)
    stop("alpha must be a finite number strictly between zero and one.")
  n <- length(p)
  if (!is.null(r.followup) &&
      (!is.numeric(r.followup) || length(r.followup) != 1L ||
       !is.finite(r.followup) || r.followup != floor(r.followup) ||
       r.followup < 1 || r.followup > n))
    stop("r.followup must be NULL or an integer between one and length(p).")
  ord <- order(inputs$log.p)
  ordered <- inputs$log.p[ord]
  raw <- vapply(seq_len(n), function(r)
    .os_combine_logp(ordered[seq.int(r, n)], settings$cutoff)$log.p.value,
    numeric(1))
  adjusted <- cummax(raw)
  tests <- data.frame(r = seq_len(n), p.value = exp(raw),
                      p.adjusted = exp(adjusted), reject = adjusted <= log(alpha),
                      log.p.value = raw, log.p.adjusted = adjusted)
  followup <- NULL
  if (!is.null(r.followup)) {
    r <- as.integer(r.followup)
    gate.reject <- raw[r] <= log(alpha)
    threshold <- if (r == n) NA_real_ else alpha / (n - r)
    log.threshold <- if (r == n) NA_real_ else log(alpha) - log(n - r)
    decisions <- if (!gate.reject) rep(FALSE, n) else if (r == n)
      rep(TRUE, n) else inputs$log.p <= log.threshold
    followup <- list(r = r, gate.p.value = exp(raw[r]),
      gate.log.p.value = raw[r], gate.reject = gate.reject,
      threshold = threshold, log.threshold = log.threshold,
      tests = data.frame(hypothesis = inputs$labels, p.value = unname(inputs$p),
                         log.p.value = unname(inputs$log.p), reject = unname(decisions)),
      plan = if (r == n)
        "Prespecified all-effects gate; no individual threshold is needed." else
        "Prespecified partial-conjunction gate followed by alpha/(n-r) tests.")
  }
  out <- list(tests = tests, lower.bound = sum(tests$reject),
    p.values = inputs$p, log.p.values = inputs$log.p,
    ordered.hypotheses = inputs$labels[ord],
    alpha = alpha, combination = method, trunc = settings$cutoff,
    count.plan = "Count-only fixed sequence; separate from the prespecified individual-follow-up family.",
    followup = followup,
    assumptions = paste("Independent true-null p-values, or the required joint",
      "stochastic bounds for every relevant intersection null; fixed combination",
      "and cutoff; r.followup must be prespecified. No simultaneous claim over",
      "an arbitrary collection of sensitivity scenarios is made."))
  class(out) <- "os_partial_conjunction"
  out
}

print.os_partial_conjunction <- function(x, digits = 4, ...) {
  dots <- list(...)
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  frame <- function(data) do.call(base::print.data.frame,
                                  c(list(x = data, digits = digits), dots))
  cat("Partial conjunction using", x$combination, "combination\n")
  frame(x$tests[, c("r", "p.value", "p.adjusted", "reject")])
  cat("Lower bound on the number of false component nulls:", x$lower.bound,
      "at alpha =", x$alpha, "\n")
  if (!is.null(x$followup)) {
    cat("Count sequence and individual follow-up are separate level-alpha plans.\n")
    cat("Prespecified follow-up target r =", x$followup$r,
        "; gate rejected:", x$followup$gate.reject, "\n")
    frame(x$followup$tests[, c("hypothesis", "p.value", "reject")])
  }
  invisible(x)
}
