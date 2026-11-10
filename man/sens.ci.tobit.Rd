\name{sens.ci.tobit}
\alias{sens.ci.tobit}
\alias{print.os_tobit_ci}
\title{Numerical inversion of matched-set Tobit-effect sensitivity tests}
\usage{
sens.ci.tobit(y, z, mset, Gamma = 1, alpha = 0.05,
  alternative = c("two-sided", "greater", "less"),
  method = c("mean", "m"), inner = 0, trim = 3, lambda = 1/2,
  grid = seq(0,1,length.out=201), tol = 1e-6)
\method{print}{os_tobit_ci}(x, ...)
}
\arguments{
\item{x}{An object returned by the corresponding function.}
\item{...}{Additional print arguments (currently unused).}
\item{y,z,mset,Gamma,method,inner,trim,lambda}{As in sens.analysis.tobit.}
\item{alpha}{Significance level strictly between zero and one.}
\item{alternative}{Direction for the reduction parameter. Greater typically gives a lower bound; less typically gives an upper bound.}
\item{grid}{Finite candidate reductions in [0,1]. Endpoints zero and one are always added. Use a finer grid to check resolution.}
\item{tol}{Positive numerical merging tolerance, no larger than .01. Roots are refined using a tighter tolerance.}
}
\description{Recomputes the censored adjusted outcomes at every candidate reduction and numerically inverts sens.analysis.tobit. It also searches for Hodges--Lehmann-type estimating roots and their sensitivity ranges.}
\details{
The search checks the grid and interval midpoints and refines detected sign
changes. Multiple accepted components and isolated accepted points are
retained. It does not assume a unique monotone crossing. Finite scans can
miss narrow components or additional crossings between checked points;
refining tol does not substitute for refining grid. These are numerical
sets based on asymptotic tests, not a certified exhaustive continuous
parameter inversion. The full permitted parameter range is [0,1], so
reported endpoints at zero or one can be model boundaries.

The estimating set contains candidates for which the statistic lies
between its minimum and maximum sensitivity expectations. At Gamma=1,
this is a zero of the centered statistic; a flat zero region is retained
rather than arbitrarily selecting one point. An all-zero response is
uninformative and gives the complete parameter range. Empty confidence or
estimating sets are returned as zero-row matrices; their envelope endpoints
are NA. A disconnected set is not replaced by its envelope in the
confidence.set or estimate.set fields. conf.int and point.estimates are
explicitly just envelopes, and may contain gaps rejected by the scan.
}
\value{An os_tobit_ci list containing confidence.set and estimate.set
(two-column component matrices), conf.int and point.estimates (their
envelopes), Gamma, alpha, alternative, method, a table of grid p-values,
grid.resolution (maximum spacing), tol, numerical=TRUE, and call.}
\examples{
y <- rep(c(.3,.5,.5),12)
z <- rep(c(1,0,0),12)
mset <- rep(1:12,each=3)
sens.ci.tobit(y,z,mset,alternative="greater",grid=seq(0,1,.05))
}
\seealso{\code{\link{sens.analysis.tobit}}, \code{\link{sens.ci.full}}}
