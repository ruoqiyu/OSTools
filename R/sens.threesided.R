sens.threesided <- function(y, z = NULL, mset = NULL, margin,
                            Gamma = 1, alpha = 0.05,
                            method = c("m", "signedrank", "alignedrank"),
                            design = c("matched", "full"), inner = 0, trim = 3,
                            lambda = 1/2, weight.par = c(1, 1, 1),
                            interval = NULL, tol = 0.01) {
  method <- match.arg(method)
  design <- match.arg(design)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(margin) || margin <= 0) stop("margin must be positive and finite, chosen before inspecting the effect.")
  if (!scalar(Gamma) || Gamma < 1) stop("Gamma must be a finite scalar at least 1.")
  if (!scalar(alpha) || alpha <= 0 || alpha >= 0.5) stop("alpha must lie between 0 and 0.5.")
  if (!scalar(tol) || tol <= 0) stop("tol must be positive and finite.")
  if (design == "matched") {
    if (!is.null(z) || !is.null(mset)) stop("For a matched design supply the response matrix; omit z and mset.")
    if ((!is.matrix(y) && !is.data.frame(y)) || nrow(y) < 2L || ncol(y) < 2L)
      stop("y must be a response matrix with at least two matched sets and two columns.")
    y <- as.matrix(y)
    if (!is.numeric(y) || any(!is.finite(y))) stop("y must contain finite numeric responses.")
    if (method == "signedrank" && ncol(y) != 2L) stop("Signed ranks require matched pairs.")
    test <- function(tau, alternative) sens.analysis(y, Gamma = Gamma, tau = tau,
        alternative = alternative, method = method, inner = inner, trim = trim,
        lambda = lambda, weight.par = weight.par)
    ci <- sens.ci(y, Gamma = Gamma, method = method, inner = inner, trim = trim,
                  lambda = lambda, weight.par = weight.par, alpha = alpha,
                  alternative = "two-sided", tol = tol, interval = interval)$Confidence.Interval
  } else {
    if (method != "m") stop("Full matching supports method='m' for three-sided inference.")
    if (is.null(z) || is.null(mset)) stop("Supply response, treatment and set-label vectors for full matching.")
    if (!is.null(interval) || !isTRUE(all.equal(weight.par, c(1, 1, 1))))
      stop("interval and nondefault weight.par are only available for a matched design.")
    test <- function(tau, alternative) sens.analysis.full(y, z, mset,
        Gamma = Gamma, method = "m", alternative = alternative, tau = tau,
        inner = inner, trim = trim, lambda = lambda)
    ci <- sens.ci.full(y, z, mset, Gamma = Gamma, inner = inner, trim = trim,
                       lambda = lambda, alpha = alpha,
                       alternative = "two-sided")$ConfidenceInterval
  }
  upper.p <- function(fit) {
    p <- if (is.list(fit$pval)) fit$pval$upperbound else fit$pval
    if (!scalar(p) || p < 0 || p > 1) stop("The analysis did not return a valid upper p-value bound.")
    p
  }
  pminus <- upper.p(test(-margin, "greater"))
  pplus <- upper.p(test(margin, "less"))
  if (!is.numeric(ci) || length(ci) != 2L || anyNA(ci) || ci[1] > ci[2])
    stop("The analysis did not return an ordered confidence interval.")
  rejected <- c(inferiority = pminus < alpha,
                equivalence = ci[2] < -margin || ci[1] > margin,
                superiority = pplus < alpha)
  tests <- data.frame(hypothesis = names(rejected),
                      boundary = c(-margin, NA_real_, margin),
                      p.value = c(pminus, NA_real_, pplus),
                      rejected = unname(rejected), row.names = NULL)
  list(tests = tests, compatible = names(rejected)[!rejected],
       equivalence.established = unname(rejected["inferiority"] && rejected["superiority"]),
       Confidence.Interval = ci, confidence.level = 1 - alpha,
       Gamma = Gamma, margin = margin, alpha = alpha, method = method, design = design)
}
