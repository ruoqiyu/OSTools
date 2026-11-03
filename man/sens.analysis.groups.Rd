\name{sens.analysis.groups}
\alias{sens.analysis.groups}
\alias{print.os_sensitivity_groups}
\title{Sensitivity Analysis in Disjoint Groups of Matched Pairs}
\usage{
sens.analysis.groups(y, group, Gamma = 1, tau = 0,
  alternative = c("two-sided", "greater", "less"),
  method = c("m", "signedrank", "mcnemar"),
  combination = c("truncated", "fisher"), trunc = 0.2, alpha = 0.05,
  inner = 0, trim = 3, lambda = 1/2, weight.par = c(1, 1, 1),
  max.groups = 15L)
\method{print}{os_sensitivity_groups}(x, digits = 4, ...)
}
\arguments{
\item{y}{Finite numeric matrix with exactly two columns: treated and control responses. Each row is an entire matched pair.}
\item{group}{One nonmissing, nonempty group label per row, aligned by position. Factors retain their level order with unused levels dropped; other labels retain first-appearance order.}
\item{Gamma}{Nonempty finite numeric vector of sensitivity bounds, each at least 1. Input order is retained.}
\item{tau}{Hypothesized constant additive effect, a finite scalar.}
\item{alternative}{Direction of the treated-minus-control difference relative to \code{tau}.}
\item{method}{Matched-pair m-test, signed-rank test, or McNemar test. McNemar requires binary outcomes and \code{tau = 0}.}
\item{combination}{Truncated-product or Fisher local intersection tests.}
\item{trunc}{Prespecified cutoff in \eqn{(0,1]}; Fisher uses 1.}
\item{alpha}{Significance level strictly between zero and one.}
\item{inner, trim, lambda, weight.par}{M-score tuning arguments passed to \code{\link{sens.analysis}}. Here \code{weight.par} must contain integers \eqn{m,m_1,m_2} with \eqn{1\leq m_1\leq m_2\leq m}.}
\item{max.groups}{Computational closure limit; see \code{\link{closed.testing}}.}
\item{x}{Object returned by \code{sens.analysis.groups}.}
\item{digits}{Printing precision.}
\item{...}{Additional arguments passed to data-frame printing.}
}
\description{Computes sensitivity upper bounds separately in disjoint groups of matched pairs, combines the group evidence, and applies closed testing to identify groups. Also reports the sensitivity bound from pooling all pairs.}
\details{
Each group's scores, scales and ranks are computed using only that group's pairs. Pooling all pairs is a separate analysis; its p-value generally differs from the combined group p-value. Groups are held fixed across the supplied sensitivity bounds. Each upper bound uses \code{\link{sens.analysis}} with the specified method and direction; direct log tails preserve small probabilities in subsequent combinations. A group with all adjusted pair differences zero has p-value 1. The m-test and signed-rank test use their existing large-sample normal approximations, which may be inadequate for small groups. McNemar bounds use binomial tails.

The assignment model requires independent assignments across matched pairs conditional on the relevant potential outcomes and covariates, with the stated within-pair sensitivity bound. Whole pairs must belong to disjoint groups. Groups may be prespecified or obtained by a discovery procedure whose predictors, unsigned response information and fitting choices are invariant to the relevant treatment-label exchanges. For a nonzero \code{tau}, invariance must hold under that tested null as well. The function neither fits a tree nor checks discovery validity. Arbitrary outcome-driven grouping does not gain validity by being supplied to this function. The method, direction and combination cutoff must satisfy the prespecification or invariance requirements for the intended analysis.

Closed testing additionally requires valid local tests for each relevant intersection; see \code{\link{closed.testing}}. Results at multiple \code{Gamma} values describe a sensitivity curve and are pointwise statements, not an additional guarantee over arbitrary searches across analyses. Rejection concerns the specified no-effect hypotheses within groups; a rejection in one group and nonrejection in another does not test whether their effects differ. Village- or hospital-level assignment is not handled by declaring individuals to be independent pairs.
}
\value{An \code{os_sensitivity_groups} list. \code{summary} has one row per supplied \code{Gamma}, with \code{pooled.p}, \code{combined.p}, corresponding rejection indicators and log probabilities. \code{tests} has one row per group per sensitivity bound, with \code{Gamma}, \code{group}, \code{pairs}, \code{raw.p}, \code{closed.p}, \code{reject}, \code{log.raw.p}, and \code{log.closed.p}. \code{closed.tests} contains the full \code{closed.testing} objects in Gamma input order. Additional fields record group labels and analysis settings and assumptions.}
\references{
Hsu, J. Y., Zubizarreta, J. R., Small, D. S. and Rosenbaum, P. R. (2015). Strong control of the familywise error rate in observational studies that discover effect modification by exploratory methods. \emph{Biometrika}, 102, 767--782. \doi{10.1093/biomet/asv034}.
}
\examples{
# Groups fixed using a baseline characteristic, independently of outcomes.
y <- cbind(c(5, 4, 6, 3, 2, 1, 4, 2), c(1, 2, 3, 1, 3, 2, 2, 4))
group <- rep(c("young", "older"), each = 4)
sens.analysis.groups(y, group, Gamma = c(1, 1.5), alternative = "greater")
}
\seealso{\code{\link{sens.analysis}}, \code{\link{closed.testing}}, \code{\link{combine.pvalues}}}
