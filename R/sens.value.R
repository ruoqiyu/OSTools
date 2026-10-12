sens.value <- function(y = NULL, D = NULL, Tobs = NULL, z = NULL, mset = NULL,
                       method = "m", design = c("matched", "full", "strat"),
                       alternative = c("greater", "less", "two-sided"),
                       alpha = 0.05, tau = 0, interval = c(1, 10),
                       tol = 1e-4, max.Gamma = 1e4, ...) {
  design <- match.arg(design)
  alternative <- match.arg(alternative)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must lie between 0 and 1.")
  if (!scalar(tau)) stop("tau must be a finite scalar.")
  if (!scalar(tol) || tol <= 0) stop("tol must be positive and finite.")
  if (!is.numeric(interval) || length(interval) != 2L || any(!is.finite(interval)) ||
      interval[1] != 1 || interval[2] <= 1)
    stop("interval must be c(1, upper), with finite upper > 1.")
  if (!scalar(max.Gamma) || max.Gamma < interval[2]) stop("max.Gamma must be finite and at least interval[2].")
  if (design != "matched" && (is.null(z) || is.null(mset))) stop("Supply z and mset for full or strat designs.")
  if (design == "matched" && (!is.null(z) || !is.null(mset))) stop("For a matched design supply the response matrix; omit z and mset.")
  if (design != "matched" && (!is.null(D) || !is.null(Tobs))) stop("D and Tobs apply only to matched McNemar analyses.")
  extra <- list(...)
  if (length(extra) && (is.null(names(extra)) || any(!nzchar(names(extra))) ||
      any(names(extra) %in% c("Gamma", "y", "D", "Tobs", "z", "mset", "stratum",
                             "method", "alternative", "tau", "level"))))
    stop("Additional arguments must be named tuning parameters, not overrides of analysis arguments.")
  evaluations <- 0L
  evaluate <- function(g) {
    evaluations <<- evaluations + 1L
    common <- list(y = y, Gamma = g, tau = tau, alternative = alternative, method = method)
    if (design == "matched") fit <- do.call(sens.analysis, c(common, list(D = D, Tobs = Tobs), extra))
    else if (design == "full") fit <- do.call(sens.analysis.full, c(common, list(z = z, mset = mset), extra))
    else fit <- do.call(sens.analysis.strat, c(common, list(z = z, stratum = mset, level = alpha), extra))
    p <- if (is.list(fit$pval)) fit$pval$upperbound else fit$pval
    if (!scalar(p) || p < 0 || p > 1) stop("The analysis did not return a finite p-value in [0,1].")
    reject <- if (design == "strat") fit$rejected else p < alpha
    if (!is.logical(reject) || length(reject) != 1L || is.na(reject)) stop("The analysis did not return a valid decision.")
    list(pval = p, rejected = reject)
  }
  result <- function(value, status, bracket, at.upper) {
    list(value = value, status = status, bracket = bracket, pval = at.upper$pval,
         rejected = at.upper$rejected, alpha = alpha, alternative = alternative,
         design = design, method = method, tau = tau, evaluations = evaluations,
         criterion = if (design == "strat") "Level-specific stratified rejection decision" else "Upper p-value bound reaches alpha")
  }
  base <- evaluate(1)
  if (!base$rejected) return(result(NA_real_, "not_significant", c(1, 1), base))
  lower <- 1
  upper <- interval[2]
  high <- evaluate(upper)
  while (high$rejected && upper < max.Gamma) {
    lower <- upper
    upper <- min(max.Gamma, 1 + 2 * (upper - 1))
    high <- evaluate(upper)
  }
  if (high$rejected) return(result(NA_real_, "above_max", c(max.Gamma, Inf), high))
  while (upper - lower > tol) {
    mid <- lower + (upper - lower)/2
    if (mid == lower || mid == upper) break
    fit <- evaluate(mid)
    if (fit$rejected) lower <- mid else {upper <- mid; high <- fit}
  }
  result(upper, "found", c(lower, upper), high)
}
