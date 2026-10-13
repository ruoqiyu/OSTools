\name{balance.check.strat}
\alias{balance.check.strat}
\title{Assess Overall and Within-Stratum Covariate Balance}
\usage{
balance.check.strat(z, X, stratum, selected = NULL, treated = 1,
  controls = NULL, denom.variance = c("pooled", "treated"),
  plot = FALSE, message = FALSE)
}
\arguments{
\item{z}{Complete group labels in the original, preselection sample.}
\item{X}{Complete covariate data frame or matrix in original row order. Numeric, logical, character and factor columns are supported. Handle missing covariates before calling this function.}
\item{stratum}{Complete exact-stratum labels, one per original observation.}
\item{selected}{Complete logical selection vector, such as \code{match.natural(...)$selected}. Default \code{NULL} assesses the full sample as both before and after.}
\item{treated}{Treated group label, default \code{1}.}
\item{controls}{Distinct control labels to assess. Default compares every other group separately with treated.}
\item{denom.variance}{Original-sample pairwise pooled variance (the average of the two group variances), or original treated variance. Uses the same denominator before and after selection.}
\item{plot}{Return the backend's ggplot balance displays as well as tables.}
\item{message}{Print the backend's balance summaries.}
}
\description{Calls the optional natstrat package to compute absolute standardized differences overall, within each exact stratum and as a control-count-weighted average across strata. Multiple controls are assessed pairwise against treated.}
\details{
Each retained assessed stratum must contain treated and control observations. Strata omitted from the design must be filtered from all inputs explicitly. Categorical covariates are expanded into indicators. For multiple controls, the original sample is restricted to the treated/control pair before standardizing, so a third group's variance cannot change that pair's diagnostic. This differs from pooling variances across all groups at once.

Standard deviations always come from the original sample. The backend uses pooled variance when the treated variance is zero. Undefined standardized differences remain undefined; do not interpret them as evidence of balance. Means alone do not establish distributional balance or eliminate hidden bias. This function does not estimate treatment effects or sampling uncertainty.
}
\value{A list with \code{comparisons}, named by control group. Each element contains the natstrat fields \code{sd_across}, \code{sd_strata}, and \code{sd_strata_avg}: tables of absolute standardized differences before and after selection. With \code{plot = TRUE}, each element also includes \code{plot_across}, \code{plot_strata}, \code{plot_strata_avg}, and \code{plot_pair}; \code{plot_strata} is a named list of plots, one per stratum. For one control, these fields are also available at the top level. Additional fields are realized \code{counts} and \code{ratios}, the \code{treated} and \code{controls} labels, and \code{denom.variance}.}
\examples{
if (requireNamespace("natstrat", quietly = TRUE)) {
  z <- rep(c(1, 0), each = 4)
  st <- rep(c("A", "B"), 4)
  X <- data.frame(x = c(1, 3, 2, 4, 0, 2, 1, 5))
  b <- balance.check.strat(z, X, st)
  b$sd_across
  b$sd_strata
}
}
\seealso{\code{\link{match.natural}}, \code{\link[natstrat]{check_balance}}}
