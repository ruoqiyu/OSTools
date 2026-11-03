sens.analysis.groups <- function(y, group, Gamma = 1, tau = 0,
    alternative = c("two-sided", "greater", "less"),
    method = c("m", "signedrank", "mcnemar"),
    combination = c("truncated", "fisher"), trunc = 0.2, alpha = 0.05,
    inner = 0, trim = 3, lambda = 1/2, weight.par = c(1, 1, 1),
    max.groups = 15L) {
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  combination <- match.arg(combination)
  settings <- .os_combination_method(combination, trunc)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!is.matrix(y) || !is.numeric(y) || ncol(y) != 2L || nrow(y) < 1L ||
      any(!is.finite(y)))
    stop("y must be a finite numeric matrix with two columns: treated, control.")
  if (!(is.character(group) || is.numeric(group) || is.logical(group) || is.factor(group)) ||
      !is.null(dim(group)) || length(group) != nrow(y) || anyNA(group))
    stop("group must supply one nonmissing group label for each pair (row of y).")
  if (is.numeric(group) && any(!is.finite(group))) stop("Numeric group labels must be finite.")
  if (anyDuplicated(as.character(unique(group))))
    stop("Distinct group labels must remain distinct when converted to character.")
  labels <- as.character(group)
  if (any(!nzchar(labels))) stop("Group labels must be nonempty.")
  groups <- if (is.factor(group)) levels(droplevels(group)) else unique(labels)
  .os_closed_settings(alpha, max.groups, length(groups))
  indices <- lapply(groups, function(g) which(labels == g))
  names(indices) <- groups
  if (!is.numeric(Gamma) || !is.null(dim(Gamma)) || !length(Gamma) ||
      any(!is.finite(Gamma)) || any(Gamma < 1))
    stop("Gamma must be a nonempty finite numeric vector with values at least one.")
  if (!scalar(tau)) stop("tau must be a finite scalar.")
  if (method == "mcnemar" && (tau != 0 || any(!y %in% c(0, 1))))
    stop("McNemar analysis requires binary responses and tau = 0.")
  if (method == "m") {
    if (!scalar(inner) || inner < 0 || !is.numeric(trim) || length(trim) != 1L ||
        is.na(trim) || trim < inner || trim <= 0 || !scalar(lambda) || lambda <= 0 || lambda >= 1)
      stop("M-score tuning requires 0 <= inner <= trim, trim > 0, and 0 < lambda < 1.")
    if (!is.numeric(weight.par) || length(weight.par) != 3L ||
        any(!is.finite(weight.par)) || any(weight.par != floor(weight.par)) ||
        weight.par[2] < 1 || weight.par[3] < weight.par[2] || weight.par[1] < weight.par[3])
      stop("weight.par must contain integers m, m1, m2 with 1 <= m1 <= m2 <= m.")
  }
  one <- function(rows, g, direction = alternative) {
    a <- y[rows, , drop = FALSE]
    d <- a[, 1] - tau - a[, 2]
    if (any(!is.finite(d))) stop("Adjusted pair differences must be finite.")
    if (all(d == 0)) return(list(pval = 1, log.pval = 0))
    if (method == "m" && direction == "two-sided") {
      hi <- one(rows, g, "greater")$log.pval
      lo <- one(rows, g, "less")$log.pval
      lp <- min(0, log(2) + min(hi, lo))
      return(list(pval = exp(lp), log.pval = lp))
    }
    fit <- sens.analysis(y = a, Gamma = g, tau = tau, alternative = direction,
      method = method, inner = inner, trim = trim, lambda = lambda, weight.par = weight.par)
    if (method == "m" && isTRUE(fit$variance == 0) &&
        isTRUE(fit$statistic == fit$expectation))
      return(list(pval = 1, log.pval = 0))
    # Calculate upper tails on the log scale before combining very small bounds.
    if (method == "m") {
      lp <- stats::pnorm(fit$deviate, lower.tail = FALSE, log.p = TRUE)
      if (direction == "two-sided") lp <- min(0, log(2) + lp)
    } else if (method == "signedrank") {
      hi <- fit$deviate$Dplus; lo <- fit$deviate$Dminus
      lp <- switch(direction,
        greater = stats::pnorm(hi, lower.tail = FALSE, log.p = TRUE),
        less = stats::pnorm(lo, log.p = TRUE),
        `two-sided` = if (hi <= 0 && lo >= 0) 0 else
          min(0, log(2) + stats::pnorm(-min(abs(hi), abs(lo)), log.p = TRUE)))
    } else {
      discordant <- sum(a[, 1] != a[, 2]); positive <- sum(a[, 1] > a[, 2])
      lp.hi <- stats::pbinom(positive - 1, discordant, g/(1 + g),
                             lower.tail = FALSE, log.p = TRUE)
      lp.lo <- stats::pbinom(positive, discordant, 1/(1 + g), log.p = TRUE)
      lp <- switch(direction, greater = lp.hi, less = lp.lo,
                   `two-sided` = min(0, log(2) + min(lp.hi, lp.lo)))
    }
    if (length(lp) != 1L || is.na(lp) || lp > 0)
      stop("The sensitivity calculation returned an invalid upper p-value bound.")
    list(pval = exp(lp), log.pval = lp)
  }
  closed <- vector("list", length(Gamma))
  pooled <- numeric(length(Gamma)); tests <- vector("list", length(Gamma))
  for (i in seq_along(Gamma)) {
    lp <- vapply(indices, function(rows) one(rows, Gamma[i])$log.pval, numeric(1))
    cl <- closed.testing(lp, method = combination, trunc = settings$cutoff,
                         alpha = alpha, log.p = TRUE, max.groups = max.groups)
    closed[[i]] <- cl
    pooled[i] <- one(seq_len(nrow(y)), Gamma[i])$log.pval
    tests[[i]] <- data.frame(Gamma = Gamma[i], group = groups, pairs = lengths(indices),
      raw.p = cl$tests$p.value, closed.p = cl$tests$p.adjusted,
      reject = cl$tests$reject, log.raw.p = cl$tests$log.p.value,
      log.closed.p = cl$tests$log.p.adjusted, row.names = NULL)
  }
  global <- vapply(closed, function(cl) cl$global.log.p.value, numeric(1))
  out <- list(summary = data.frame(Gamma = Gamma, pooled.p = exp(pooled),
      combined.p = exp(global), pooled.reject = pooled <= log(alpha),
      combined.reject = global <= log(alpha), log.pooled.p = pooled,
      log.combined.p = global),
    tests = do.call(rbind, tests), closed.tests = closed, group = labels,
    Gamma = Gamma, tau = tau, alternative = alternative, method = method,
    combination = combination, trunc = settings$cutoff, alpha = alpha,
    assumptions = paste("Independent assignments across matched pairs under the stated Gamma bound;",
      "groups fixed in advance or obtained using an appropriate assignment-invariant discovery procedure;",
      "valid local tests for every relevant intersection. m and signedrank use normal approximations.",
      "Results are pointwise in Gamma and are no-effect tests, not tests of equal subgroup effects."))
  class(out) <- "os_sensitivity_groups"
  out
}

print.os_sensitivity_groups <- function(x, digits = 4, ...) {
  dots <- list(...)
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  frame <- function(z) do.call(base::print.data.frame, c(list(x = z, digits = digits), dots))
  cat("Matched-pair grouped sensitivity analysis using", x$method, "tests\n")
  frame(x$summary[, c("Gamma", "pooled.p", "combined.p", "combined.reject")])
  cat("Closed-testing subgroup decisions:\n")
  frame(x$tests[, c("Gamma", "group", "pairs", "raw.p", "closed.p", "reject")])
  invisible(x)
}
