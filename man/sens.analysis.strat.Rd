\name{sens.analysis.strat}
\alias{sens.analysis.strat}
\title{Sensitivity Analysis for General Continuous-Outcome Strata}
\usage{
sens.analysis.strat(y, z, stratum, Gamma = 1, tau = 0,
  alternative = c("two-sided", "greater", "less"),
  method = c("m", "alignedrank", "scores"), inner = 0, trim = 3,
  lambda = 1/2, align = c("mean", "median", "huber", "hl"),
  level = 0.05, bound.method = c("BU", "RK", "LS", "AD"),
  detail = FALSE)
}
\arguments{
\item{y}{Finite numeric response vector, or precomputed scores for method = "scores".}
\item{z}{One 0-1 treatment indicator per response.}
\item{stratum}{One nonmissing stratum label per response; character, numeric and factor labels are accepted.}
\item{Gamma}{Finite sensitivity parameter at least 1, not its logarithm.}
\item{tau}{Hypothesized constant additive treatment effect. Adjusted responses are y - tau * z. Must be zero for precomputed scores.}
\item{alternative}{Two-sided, increasing-response, or decreasing-response alternative. Default is two-sided.}
\item{method}{Huber M-scores, pooled aligned ranks, or supplied scores.}
\item{inner}{Nonnegative inner trimming threshold for M-scores.}
\item{trim}{Outer trimming threshold greater than inner; Inf gives untrimmed scores.}
\item{lambda}{Scale quantile strictly between zero and one for M-scores.}
\item{align}{Location used within strata before ranking: mean, median, Huber estimate, or Hodges--Lehmann estimate. Used only for aligned ranks.}
\item{level}{Significance level strictly between zero and one.}
\item{bound.method}{Moment calculation passed to \code{senstrat::senstrat}: BU, RK, LS or AD. LS requires large samples within strata. Default BU uses exact moments in the Normal approximation.}
\item{detail}{Whether to retain complete directional backend results.}
}
\description{Wraps the senstrat package for strata with arbitrary numbers of treated and control subjects, including propensity-score subclasses.}
\details{Strata containing only one treatment group are excluded before scoring. Scores and the scale therefore depend only on retained strata. M-scores and aligned ranks are computed using senstrat's score functions. Supplied scores are left unchanged. Labels are recoded internally to consecutive integers. A single retained stratum is allowed, with the backend's recommendation to use its two-sample routine emitted as a warning.

Directional decisions use the linear Taylor criterion at the specified level. Two-sided decisions combine two directional tests, each at level/2. The approximate two-sided p-value doubles the smaller directional approximate p-value, capped at one. Directional reported Normal probabilities are evaluated directly from the backend deviates to retain small, representable tails. The p-value is descriptive; the level-specific rejection decision is authoritative and is what \code{sens.value} uses. The calculation uses a Normal approximation, not a finite-sample exact test. Constant scores within every retained stratum return p-value 1 and no rejection.

Stratified scores and weighting differ from some full-matching M-statistics. This function does not replace \code{sens.analysis.full} for its matched-set design.}
\value{A list containing \code{pval} (approximate p-value), \code{rejected} (decision at the stated level), \code{decision.margin} (linear Taylor criterion; nonpositive means rejection for nondegenerate data), \code{Gamma}, \code{tau}, \code{alternative}, \code{level}, \code{method}, \code{bound.method}, \code{directional} (direction-specific results), \code{scores} (in original row order, with NA for excluded rows), \code{excluded} (excluded row indices), and \code{pval.type}. For completely constant within-stratum scores, decision.margin is set to 1 to indicate no rejection.}
\examples{
y <- c(8, 7, 4, 3, 9, 6, 5, 2, 10, 8, 6, 4)
z <- rep(c(1, 1, 0, 0), 3)
s <- rep(c("a", "b", "c"), each = 4)
sens.analysis.strat(y, z, s, Gamma = 2, alternative = "greater")
sens.analysis.strat(y, z, s, Gamma = 2, method = "alignedrank")
}
\references{Rosenbaum, P. R. (2018). Sensitivity analysis for stratified comparisons in an observational study of the effect of smoking on homocysteine levels. The Annals of Applied Statistics, 12, 2312--2334. \doi{10.1214/18-AOAS1153}.}
\seealso{\code{\link{sens.value}}, \code{\link{sens.analysis.full}}, \code{senstrat::senstrat}}
