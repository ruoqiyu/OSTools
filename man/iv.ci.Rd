\name{iv.ci}
\alias{iv.ci}
\alias{print.os_iv_ci}
\title{Grid Confidence Sets for Rank-Based Instrumental-Variable Inference}
\description{Inverts two-sided adjusted-outcome rank tests over a finite grid, preserving all accepted components rather than replacing rejected gaps by an interval.}
\usage{
iv.ci(y, treatment=NULL, instrument=NULL, mset=NULL,
  method=c("signedrank", "ranksum"),
  Gamma=1, level=0.95, grid=NULL, interval=c(-10,10), tol=0.01,
  exact=NULL, correct=TRUE, zero.method=c("wilcox", "pratt"),
  d=NULL, z=NULL)
\method{print}{os_iv_ci}(x, ...)
}
\arguments{
\item{y, treatment, instrument, mset, method, Gamma, exact, correct, zero.method}{As in \code{\link{iv.test}}. Matrix inputs describe pairs directly; vector inputs need two instrument levels for rank sums or pair labels for signed ranks.}
\item{d, z}{Legacy aliases: d is treatment, z is instrument. Their meanings do not change. Conflicting canonical and legacy names are rejected.}
\item{level}{One confidence level strictly between 0 and 1.}
\item{grid}{Optional vector of at least two strictly increasing finite candidate effects. When supplied, it replaces \code{interval} and \code{tol}.}
\item{interval}{Two increasing finite search limits used when \code{grid=NULL}.}
\item{tol}{Positive grid spacing used when \code{grid=NULL}; it is a spacing, not a root-finding tolerance.}
\item{x}{A returned confidence-set object.}
\item{...}{Additional arguments to print the components.}
}
\details{
In Chapter 15 notation, supply \code{y=Y, treatment=Z, instrument=A}. The positional argument order is unchanged.

An effect is accepted when its two-sided p-value, or upper sensitivity p-value bound, is at least one minus the confidence level (up to relative floating-point tolerance). All consecutive runs of accepted grid points are returned separately. Their endpoints are grid approximations, not interpolated confidence limits. A narrow component or rejected gap between grid points can be missed; use suitable spacing and check stability.

Boundary contact warns that the search should be extended. It does not prove that the full confidence set is unbounded, and a reported boundary value must not be treated as an established finite confidence limit. No accepted grid points also triggers a warning; this does not establish that no effects would be accepted outside the supplied grid. Weak-instrument confidence sets may be unbounded, empty or disconnected, so these distinctions matter.

The reported estimate is a grid value maximizing the Gamma-1 two-sided p-value. It is an approximate rank-based point estimate, not a Wald ratio. The first maximizer is returned for a plateau; if all Gamma-1 p-values are equal, the point estimate is reported as missing. Sensitivity sets retain this Gamma-1 estimate. Causal interpretations and rank-test assumptions are described in \code{\link{iv.test}}.
}
\value{An \code{os_iv_ci} list with:
\item{estimate}{Gamma-1 grid p-value maximizer, or NA for a flat p-value curve.}
\item{components}{Numeric two-column matrix of accepted grid components, with columns lower and upper; zero rows if none were accepted.}
\item{span}{Smallest and largest accepted grid values, or two NAs when empty. A span may contain rejected effects.}
\item{touches.boundary, boundary.contact}{Logical boundary warning flag and named lower/upper contact indicators.}
\item{empty, disconnected}{Whether the accepted grid set is empty or has multiple components.}
\item{grid, p.values, accepted}{Candidate effects, upper p-values, and acceptance indicators.}
\item{level, Gamma}{Confidence level and sensitivity parameter.}
\item{first.stage, n, n.discarded}{Observed treatment difference, group counts and unmatched-subject count.}
\item{method, exact, correct, zero.method, call}{Analysis settings and call; exact records the requested setting.}
}
\seealso{\code{\link{iv.test}}, \code{\link{sens.ci.iv}}}
\examples{
y <- cbind(c(2,3,4,5,6,7),rep(0,6))
Z <- cbind(rep(1,6),rep(0,6))
fit <- iv.ci(y,treatment=Z,interval=c(-10,15),tol=0.1,exact=FALSE)
fit
fit$components
}
