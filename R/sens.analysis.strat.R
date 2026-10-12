sens.analysis.strat <- function(y, z, stratum, Gamma = 1, tau = 0,
                                alternative = c("two-sided", "greater", "less"),
                                method = c("m", "alignedrank", "scores"),
                                inner = 0, trim = 3, lambda = 1/2,
                                align = c("mean", "median", "huber", "hl"),
                                level = 0.05, bound.method = c("BU", "RK", "LS", "AD"),
                                detail = FALSE) {
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  align <- match.arg(align)
  bound.method <- match.arg(bound.method)
  if (!is.numeric(y) || !is.null(dim(y)) || !length(y) || any(!is.finite(y)))
    stop("y must be a finite numeric response vector.")
  if (length(z) != length(y) || anyNA(z) || any(!z %in% c(0, 1)))
    stop("z must contain one 0-1 treatment indicator per response.")
  if (length(stratum) != length(y) || anyNA(stratum))
    stop("stratum must contain one nonmissing label per response.")
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(Gamma) || Gamma < 1) stop("Gamma must be a finite scalar at least 1.")
  if (!scalar(tau)) stop("tau must be a finite scalar.")
  if (!scalar(level) || level <= 0 || level >= 1) stop("level must lie between 0 and 1.")
  if (!is.logical(detail) || length(detail) != 1L || is.na(detail)) stop("detail must be TRUE or FALSE.")
  if (method == "scores" && tau != 0) stop("For method='scores', adjust and score responses first and set tau=0.")
  if (method == "m" && (!scalar(inner) || inner < 0 || !is.numeric(trim) ||
      length(trim) != 1L || is.na(trim) || trim <= inner ||
      !scalar(lambda) || lambda <= 0 || lambda >= 1))
    stop("M-score tuning requires 0 <= inner < trim and 0 < lambda < 1.")
  st <- as.integer(factor(stratum))
  groups <- split(seq_along(y), st)
  mixed <- vapply(groups, function(i) any(z[i] == 1) && any(z[i] == 0), logical(1))
  keep <- st %in% as.integer(names(groups)[mixed])
  if (!any(keep)) stop("At least one stratum must contain both treatment groups.")
  excluded <- which(!keep)
  original.y <- y
  y <- y[keep]; z <- z[keep]; st <- as.integer(factor(st[keep]))
  adjusted <- y - tau * z
  constant <- all(vapply(split(adjusted, st), function(v) max(v) == min(v), logical(1)))
  if (method == "m" && constant) scores <- rep(0, length(y))
  else if (method == "m") scores <- senstrat::mscores(y, z, st = st,
      inner = inner, trim = trim, lambda = lambda, tau = tau)
  else if (method == "alignedrank") scores <- senstrat::hodgeslehmann(y, z, st,
      align = align, tau = tau)
  else scores <- y
  if (any(!is.finite(scores))) stop("Scoring produced nonfinite values; check the scale and tuning parameters.")
  degenerate <- all(vapply(split(scores, st), function(v) max(v) == min(v), logical(1)))
  one <- function(direction, test.level) {
    if (degenerate) return(list(pval = 1, rejected = FALSE, decision.margin = 1,
                                result = NULL, backend = NULL))
    fit <- senstrat::senstrat(scores, z, st, gamma = Gamma,
                             alternative = direction, level = test.level,
                             method = bound.method, detail = TRUE)
    # Preserve backend decisions, but evaluate reported tails directly.
    upper.tail <- direction == "greater"
    fit$LinearBoundResult["P-value"] <- stats::pnorm(
      fit$LinearBoundResult["Deviate"], lower.tail = !upper.tail)
    fit$Separable["P-value"] <- stats::pnorm(
      fit$Separable["Deviate"], lower.tail = !upper.tail)
    margin <- unname(fit$lambda["Linear Taylor bound"])
    list(pval = unname(fit$LinearBoundResult["P-value"]),
         rejected = margin <= 0, decision.margin = margin,
         result = fit$LinearBoundResult, backend = if (detail) fit else NULL)
  }
  if (alternative == "two-sided") {
    directional <- list(greater = one("greater", level/2), less = one("less", level/2))
    pval <- min(1, 2 * min(directional$greater$pval, directional$less$pval))
    rejected <- directional$greater$rejected || directional$less$rejected
    decision.margin <- min(directional$greater$decision.margin, directional$less$decision.margin)
  } else {
    directional <- stats::setNames(list(one(alternative, level)), alternative)
    pval <- directional[[1]]$pval
    rejected <- directional[[1]]$rejected
    decision.margin <- directional[[1]]$decision.margin
  }
  output.scores <- rep(NA_real_, length(original.y)); output.scores[keep] <- scores
  list(pval = pval, rejected = rejected, decision.margin = decision.margin,
       Gamma = Gamma, tau = tau, alternative = alternative, level = level,
       method = method, bound.method = bound.method,
       directional = directional, scores = output.scores, excluded = excluded,
       pval.type = "Approximate; use rejected for the decision at the stated level.")
}
