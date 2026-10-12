sens.amplify <- function(Gamma, Lambda = NULL) {
  if (!is.numeric(Gamma) || !length(Gamma) || any(!is.finite(Gamma)) || any(Gamma < 1))
    stop("Gamma must contain finite numbers at least 1.")
  symmetric <- Gamma + sqrt(Gamma - 1) * sqrt(Gamma + 1)
  if (is.null(Lambda)) {
    Lambda <- symmetric
    Delta <- symmetric
  } else {
    if (!is.numeric(Lambda) || !length(Lambda) || any(!is.finite(Lambda)))
      stop("Lambda must contain finite numbers.")
    n <- max(length(Gamma), length(Lambda))
    if (!length(Gamma) %in% c(1L, n) || !length(Lambda) %in% c(1L, n))
      stop("Gamma and Lambda must have equal lengths, or one must be a scalar.")
    Gamma <- rep(Gamma, length.out = n)
    Lambda <- rep(Lambda, length.out = n)
    symmetric <- rep(symmetric, length.out = n)
    if (any(Lambda < 1) || any(Gamma > 1 & Lambda <= Gamma))
      stop("Lambda must exceed Gamma when Gamma > 1, and be at least 1 when Gamma = 1.")
    Delta <- rep(1, n)
    use <- Gamma > 1
    Delta[use] <- Gamma[use] + (Gamma[use] - 1) *
      ((Gamma[use] + 1)/(Lambda[use] - Gamma[use]))
  }
  data.frame(Gamma = Gamma, Lambda = Lambda, Delta = Delta,
             probability.lower = 1/(1 + Gamma),
             probability.upper = 1 - 1/(1 + Gamma),
             symmetric = symmetric)
}
