\name{sens.threesided}
\alias{sens.threesided}
\title{Three-Sided Sensitivity Inference with an Equivalence Margin}
\usage{
sens.threesided(y, z = NULL, mset = NULL, margin, Gamma = 1,
  alpha = 0.05, method = c("m", "signedrank", "alignedrank"),
  design = c("matched", "full"), inner = 0, trim = 3,
  lambda = 1/2, weight.par = c(1, 1, 1), interval = NULL, tol = 0.01)
}
\arguments{
\item{y}{Matched response matrix, with treated responses in column one, or full-matching response vector. Matched matrices need at least two rows.}
\item{z}{Treatment vector for full matching; omit for matched matrices.}
\item{mset}{Set-label vector for full matching; omit for matched matrices.}
\item{margin}{Positive equivalence margin on the response scale, specified before inspecting the effect estimate.}
\item{Gamma}{Finite sensitivity parameter at least 1.}
\item{alpha}{Significance level strictly between zero and 0.5.}
\item{method}{M-statistic, signed rank, or aligned rank for matched designs. Signed ranks require pairs. Full matching supports M-statistics.}
\item{design}{Matched matrix or full-matching vectors.}
\item{inner}{M-score inner trimming parameter.}
\item{trim}{M-score outer trimming parameter.}
\item{lambda}{M-score scale quantile.}
\item{weight.par}{M-statistic weight parameters for matched matrices; full matching requires the default.}
\item{interval}{Optional confidence-interval search range for matched matrices only.}
\item{tol}{Matched confidence-interval root tolerance, or signed-rank grid spacing, default 0.01. Smaller values give finer numerical resolution. Full matching uses the backend's interval algorithm.}
}
\description{Tests inferiority, equivalence and superiority regions using the sensitivity analyses already provided by OSTools.}
\details{Under the constant additive effect model, inferiority is tau < -margin, equivalence is -margin <= tau <= margin, and superiority is tau > margin. Inferiority is tested at the boundary -margin against larger effects, and superiority at margin against smaller effects, each at level alpha. The equivalence null is rejected only if the two-sided (1-alpha) sensitivity interval has no overlap with the closed equivalence region.

Evidence for equivalence requires rejection of both inferiority and superiority. Failure to reject the equivalence null alone is not evidence for equivalence. The two one-sided tests correspond to a (1-2*alpha) interval being inside the margin; the returned (1-alpha) interval serves the separate test of the equivalence null. The method inherits the approximations and additive-effect assumptions of its underlying tests and intervals. Signed-rank calculations stop if the search range may truncate the interval.}
\value{A list with \code{tests}, a three-row data frame containing hypothesis, boundary, p.value and rejected; \code{compatible}, the effect regions not rejected; \code{equivalence.established}; \code{Confidence.Interval}; \code{confidence.level}; \code{Gamma}; \code{margin}; \code{alpha}; \code{method}; and \code{design}. The equivalence-null p.value is NA because its decision comes from interval inversion.}
\examples{
difference <- seq(-0.3, 0.3, length.out = 30)
control <- seq(0, 2, length.out = 30)
y <- cbind(control + difference, control)
sens.threesided(y, margin = 0.5)
sens.threesided(y, margin = 0.5, Gamma = 2)
}
\seealso{\code{\link{sens.analysis}}, \code{\link{sens.ci}}, \code{\link{sens.analysis.full}}, \code{\link{sens.ci.full}}}
