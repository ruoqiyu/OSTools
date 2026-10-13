contrast.test <- function(y = NULL, group = NULL, contrast, means = NULL,
                          se = NULL, n = NULL, null = 0,
                          alternative = c("two-sided", "greater", "less"),
                          conf.level = 0.95, na.rm = FALSE) {
  if (identical(alternative, "two.sided")) alternative <- "two-sided"
  alternative <- match.arg(alternative)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(null)) stop("null must be a finite numeric scalar.")
  if (!scalar(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be strictly between zero and one.")
  groups <- .contrast.groups(y, group, means, se, n, na.rm)
  a <- .contrast.align(contrast, groups$labels, "contrast")
  if (!any(a != 0)) stop("contrast must contain at least one nonzero coefficient.")
  estimate <- sum(a * groups$means)
  delta <- estimate - null
  if (!is.finite(estimate) || !is.finite(delta))
    stop("The contrast estimate overflowed; rescale the responses or coefficients.")

  # Scale before squaring: tiny/large standard errors should not destroy df.
  contributions <- abs(a) * groups$se
  scale <- max(contributions)
  if (!is.finite(scale) || scale <= 0)
    stop("The contrast must have a positive finite standard error; rescale if needed.")
  v <- (contributions / scale)^2
  total <- sum(v)
  stderr <- scale * sqrt(total)
  df <- total^2 / sum(v^2 / (groups$n - 1))
  if (!is.finite(stderr) || stderr <= 0 || !is.finite(df) || df <= 0)
    stop("The standard error or degrees of freedom is not finite and positive.")
  statistic <- delta / stderr
  # Direct log tails retain small p-values and describe probability underflow.
  log.p <- switch(alternative,
    greater = stats::pt(statistic, df, lower.tail = FALSE, log.p = TRUE),
    less = stats::pt(statistic, df, lower.tail = TRUE, log.p = TRUE),
    `two-sided` = min(0, log(2) + stats::pt(abs(statistic), df,
                                          lower.tail = FALSE, log.p = TRUE)))
  tail <- if (alternative == "two-sided") (1 - conf.level) / 2 else 1 - conf.level
  width <- stats::qt(tail, df, lower.tail = FALSE) * stderr
  ci <- switch(alternative,
    greater = c(estimate - width, Inf),
    less = c(-Inf, estimate + width),
    `two-sided` = c(estimate - width, estimate + width))
  required <- switch(alternative, greater = 1L, less = 2L, `two-sided` = 1:2)
  if (any(!is.finite(ci[required])) || anyNA(ci) || ci[1] > ci[2])
    stop("The confidence interval overflowed or is invalid; rescale the responses.")
  attr(ci, "conf.level") <- conf.level
  result <- list(statistic = c(t = statistic), parameter = c(df = df),
                 p.value = exp(log.p), log.p.value = log.p,
                 conf.int = ci, estimate = c(contrast = estimate),
                 null.value = c(contrast = null), stderr = stderr,
                 alternative = if (alternative == "two-sided") "two.sided" else alternative,
                 method = "Welch test for an independent-group linear contrast",
                 data.name = if (groups$raw) "responses grouped by group" else "independent-group summary statistics",
                 coefficients = a, group.summary = groups$summary)
  class(result) <- "htest"
  result
}
