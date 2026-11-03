combine.pvalues <- function(p, method = c("fisher", "truncated"),
                            trunc = 0.2, log.p = FALSE) {
  method <- match.arg(method)
  settings <- .os_combination_method(method, trunc)
  inputs <- .os_pvalues(p, log.p)
  result <- .os_combine_logp(inputs$log.p, settings$cutoff)
  fisher <- settings$cutoff == 1
  out <- list(
    statistic = if (fisher) c(`X-squared` = -2 * result$log.product) else
      c(`log(product)` = result$log.product),
    parameter = if (fisher) c(df = 2 * length(p)) else
      c(n = length(p), trunc = settings$cutoff),
    p.value = exp(result$log.p.value),
    log.p.value = result$log.p.value,
    alternative = "at least one component null hypothesis is false",
    method = if (fisher) "Fisher combination of p-values" else
      "Truncated-product combination of p-values",
    data.name = paste(deparse(substitute(p)), collapse = " "),
    p.values = inputs$p, log.p.values = inputs$log.p,
    combination = method, trunc = settings$cutoff,
    n.retained = result$retained,
    assumptions = paste("Independent null p-values, or joint stochastic bounds",
      "valid for this monotone combination. These assumptions are not checked from p."))
  class(out) <- "htest"
  out
}
