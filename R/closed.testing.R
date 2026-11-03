closed.testing <- function(p, method = c("truncated", "fisher"), trunc = 0.2,
                           alpha = 0.05, log.p = FALSE, max.groups = 15L) {
  method <- match.arg(method)
  settings <- .os_combination_method(method, trunc)
  inputs <- .os_pvalues(p, log.p)
  .os_closed_settings(alpha, max.groups, length(p))
  n <- length(p)
  subsets <- unlist(lapply(seq_len(n), function(k)
    utils::combn(seq_len(n), k, simplify = FALSE)), recursive = FALSE)
  # Singleton tests use their original p-values, even above the truncation cutoff.
  local <- vapply(subsets, function(k) {
    if (length(k) == 1L) inputs$log.p[k] else
      .os_combine_logp(inputs$log.p[k], settings$cutoff)$log.p.value
  }, numeric(1))
  adjusted <- vapply(seq_len(n), function(j)
    max(local[vapply(subsets, function(k) j %in% k, logical(1))]), numeric(1))
  tests <- data.frame(hypothesis = inputs$labels, p.value = unname(inputs$p),
    p.adjusted = exp(adjusted), reject = adjusted <= log(alpha),
    log.p.value = unname(inputs$log.p), log.p.adjusted = adjusted)
  intersections <- data.frame(id = seq_along(subsets),
    hypothesis = vapply(subsets, function(k)
      paste(inputs$labels[k], collapse = " & "), character(1)),
    size = lengths(subsets), p.value = exp(local), log.p.value = local,
    local.reject = local <= log(alpha))
  intersections$indices <- I(subsets)
  out <- list(tests = tests, intersections = intersections,
    global.p.value = exp(local[length(local)]), global.log.p.value = local[length(local)],
    alpha = alpha, combination = method, trunc = settings$cutoff,
    assumptions = paste("Every intersection requires a valid local test:",
      "independent true-null p-values, or the required joint stochastic bounds.",
      "The combination and cutoff must be fixed independently of the signed results.",
      "The function does not check these conditions or validate a discovery procedure."))
  class(out) <- "os_closed_testing"
  out
}

.os_closed_settings <- function(alpha, max.groups, n) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1)
    stop("alpha must be a finite number strictly between zero and one.")
  if (!is.numeric(max.groups) || length(max.groups) != 1L ||
      !is.finite(max.groups) || max.groups != floor(max.groups) ||
      max.groups < 1 || max.groups > 20)
    stop("max.groups must be an integer from one to 20.")
  if (n > max.groups)
    stop("The number of groups exceeds max.groups; closure requires 2^G - 1 local tests.")
  invisible(NULL)
}

print.os_closed_testing <- function(x, digits = 4, ...) {
  dots <- list(...)
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  frame <- function(z) do.call(base::print.data.frame, c(list(x = z, digits = digits), dots))
  cat("Closed testing using", x$combination, "combination\n")
  frame(x$tests[, c("hypothesis", "p.value", "p.adjusted", "reject")])
  cat("Joint-null p-value:", format(x$global.p.value, digits = digits), "\n")
  invisible(x)
}
