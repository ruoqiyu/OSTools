balance.check.strat <- function(z, X, stratum, selected = NULL, treated = 1,
                                controls = NULL, denom.variance = c("pooled", "treated"),
                                plot = FALSE, message = FALSE) {
  .natural.require()
  if (is.matrix(X)) X <- as.data.frame(X)
  input <- .natural.inputs(z, stratum, X, treated)
  z <- input$z; stratum <- input$stratum; treated <- input$treated
  denom.variance <- match.arg(denom.variance)
  .natural.flag(plot, "plot"); .natural.flag(message, "message")
  if (is.null(selected)) selected <- rep(TRUE, nrow(X))
  if (!is.logical(selected) || length(selected) != nrow(X) || anyNA(selected))
    stop("selected must be a complete logical vector with one entry per row.")
  if (is.null(controls)) controls <- levels(z)[levels(z) != treated]
  controls <- as.character(controls)
  if (!length(controls) || anyNA(controls) || anyDuplicated(controls) ||
      any(!controls %in% levels(z)) || treated %in% controls)
    stop("controls must be distinct observed groups other than treated.")
  for (i in seq_along(X)) {
    if (anyNA(X[[i]])) stop("X must be complete; handle missing covariates before assessing balance.")
    if (is.character(X[[i]])) X[[i]] <- factor(X[[i]])
    if (!is.numeric(X[[i]]) && !is.logical(X[[i]]) && !is.factor(X[[i]]))
      stop("X columns must be numeric, logical, character or factor covariates.")
    if (is.numeric(X[[i]]) && any(is.infinite(X[[i]]))) stop("Numeric covariates cannot contain infinite values.")
  }
  assessed <- z %in% c(treated, controls)
  counts <- table(droplevels(z[assessed & selected]), stratum[assessed & selected])
  if (!all(c(treated, controls) %in% rownames(counts)) || any(counts < 1))
    stop("Each assessed stratum must retain treated and control observations; filter omitted strata first.")
  comparisons <- lapply(controls, function(control) {
    pair <- z %in% c(treated, control)
    covariates <- X[pair, , drop = FALSE]
    # natstrat 2.0.1 simplifies a one-covariate variance matrix to a vector.
    # A duplicate column preserves its dimensions; remove it from all returns.
    duplicate <- NULL
    if (ncol(covariates) == 1L && !is.factor(covariates[[1]])) {
      duplicate <- make.unique(c(names(covariates), "OSToolsDuplicateCovariate"))[2]
      covariates[[duplicate]] <- covariates[[1]]
    }
    out <- .natural.backend(natstrat::check_balance(z = droplevels(z[pair]), X = covariates,
      st = stratum[pair], selected = selected[pair], treated = treated,
      control = control, denom_variance = denom.variance, plot = plot, message = message))
    if (!is.null(duplicate)) {
      for (field in c("sd_across", "sd_strata_avg"))
        out[[field]] <- out[[field]][rownames(out[[field]]) != duplicate, , drop = FALSE]
      out$sd_strata <- out$sd_strata[out$sd_strata$covariate != duplicate, , drop = FALSE]
      strip <- function(p) {
        if (inherits(p, "ggplot")) {
          p$data <- droplevels(p$data[as.character(p$data$covariates) != duplicate, , drop = FALSE])
          p
        } else lapply(p, strip)
      }
      if (plot) for (field in c("plot_across", "plot_strata", "plot_strata_avg", "plot_pair"))
        out[[field]] <- strip(out[[field]])
    }
    out
  })
  names(comparisons) <- controls
  result <- list(comparisons = comparisons, counts = counts,
                 ratios = sweep(counts, 2, counts[treated, ], "/"),
                 treated = treated, controls = controls, denom.variance = denom.variance)
  # Binary callers can use the backend's existing balance/plot field names.
  if (length(controls) == 1L) result <- c(comparisons[[1]], result)
  result
}
