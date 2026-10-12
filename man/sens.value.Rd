\name{sens.value}
\alias{sens.value}
\title{Find the Sensitivity Boundary for a Specified Test}
\usage{
sens.value(y = NULL, D = NULL, Tobs = NULL, z = NULL, mset = NULL,
  method = "m", design = c("matched", "full", "strat"),
  alternative = c("greater", "less", "two-sided"), alpha = 0.05,
  tau = 0, interval = c(1, 10), tol = 1e-4, max.Gamma = 1e4, ...)
}
\arguments{
\item{y}{Response matrix for matched designs, or response vector for full or strat designs.}
\item{D}{Number of discordant pairs for the matched McNemar test.}
\item{Tobs}{Observed number of positive discordant pairs for the matched McNemar test.}
\item{z}{Treatment vector, required for full and strat designs.}
\item{mset}{Matched-set or stratum labels, required for full and strat designs.}
\item{method}{Method supported by the selected design's analysis function.}
\item{design}{Use \code{sens.analysis}, \code{sens.analysis.full}, or \code{sens.analysis.strat}.}
\item{alternative}{Specified alternative; the default is greater. The sensitivity value depends on this choice.}
\item{alpha}{Significance level strictly between zero and one.}
\item{tau}{Hypothesized additive effect, default zero.}
\item{interval}{Initial search bracket c(1, upper), with upper greater than 1. Expanded automatically if necessary.}
\item{tol}{Positive absolute tolerance for the final Gamma bracket.}
\item{max.Gamma}{Finite maximum permitted Gamma, at least the initial upper limit.}
\item{...}{Named tuning arguments passed to the analysis function, for example inner, trim, lambda, weight.par, align or bound.method. Analysis arguments cannot be overridden.}
}
\description{Finds the transition from rejection to nonrejection as the allowed hidden bias increases.}
\details{Bisection uses rejection decisions, allowing nonsmooth score calculations. Matched and full designs reject when their upper p-value bound is below alpha. Stratified designs use the level-specific rejection decision rather than treating an approximate p-value as an exact decision rule. The search starts at Gamma = 1. It doubles the distance of the upper bracket from 1 until nonrejection is found or max.Gamma is reached. The reported value is the nonrejecting end of the final bracket, and its accuracy is limited by the underlying test approximation and tol. This search assumes a monotone rejection-to-nonrejection transition in Gamma, as intended by these sensitivity analyses. No sensitivity value is claimed when the original test does not reject or the transition has not been bracketed.}
\value{A list with \code{value} (the boundary estimate, or NA), \code{status} ("found", "not_significant", or "above_max"), \code{bracket}, \code{pval} and \code{rejected} at the reported upper limit, \code{alpha}, \code{alternative}, \code{design}, \code{method}, \code{tau}, \code{evaluations}, and \code{criterion}. For above_max, bracket is c(max.Gamma, Inf), indicating a lower bound rather than an infinite sensitivity value. For not_significant, bracket is c(1,1).}
\examples{
sens.value(D = 122, Tobs = 110, method = "mcnemar",
           alternative = "greater")
sens.value(D = 122, Tobs = 110, method = "mcnemar",
           alternative = "two-sided")
# A result with status above_max reports censoring, not a fabricated root.
sens.value(D = 122, Tobs = 110, method = "mcnemar",
           interval = c(1, 2), max.Gamma = 2)
}
\seealso{\code{\link{sens.analysis}}, \code{\link{sens.analysis.full}}, \code{\link{sens.analysis.strat}}, \code{\link{sens.amplify}}}
