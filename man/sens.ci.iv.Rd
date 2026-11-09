\name{sens.ci.iv}
\alias{sens.ci.iv}
\title{Sensitivity Confidence Sets for Paired Instrumental Variables}
\description{A backward-compatible paired interface to the grid inversion in \code{iv.ci}.}
\usage{
sens.ci.iv(encouraged.treatment, unencouraged.treatment,
  encouraged.response, unencouraged.response, Gamma=1, alpha=0.05,
  tol=0.01, interval=c(-10,10), details=FALSE)
}
\arguments{
\item{encouraged.treatment, unencouraged.treatment}{Finite treatment-received vectors of identical nonzero length, ordered by pair.}
\item{encouraged.response, unencouraged.response}{Corresponding finite response vectors of that same length.}
\item{Gamma}{Finite sensitivity parameter at least 1.}
\item{alpha}{One minus the confidence level, strictly between 0 and 1.}
\item{tol}{Positive grid spacing.}
\item{interval}{Increasing finite pair of search-grid endpoints.}
\item{details}{If TRUE, return every accepted component and diagnostics as an \code{os_iv_ci} object. The default preserves the historical two-endpoint interface.}
}
\details{
Chapter 15 denotes the instrument by A and treatment received by Z. The encouraged and unencouraged treatment arguments contain Z, ordered by the instrument A. This positional interface is unchanged.

This interface uses normal signed-rank approximations with no continuity correction and Pratt zero handling, matching the historical moment calculation. It delegates to \code{iv.ci}, which now checks paired lengths, finite data, sensitivity/confidence settings and the grid. The pair can be oriented either way; inversion no longer obtains separate endpoints under a presumed first-stage sign.

Boundary contact and an empty accepted grid issue warnings. With disconnected sets, the legacy endpoint vector spans rejected gaps and issues a warning; use \code{details=TRUE} to retain the actual components. These are finite-grid calculations, which cannot by themselves establish globally empty or unbounded confidence sets. See \code{\link{iv.ci}} for limitations and \code{\link{iv.test}} for the identifying assumptions. Gamma bounds hidden bias in instrument assignment and does not relax exclusion.
}
\value{By default, an un-named numeric vector of the first and last accepted grid values (two NAs if no grid values are accepted). With details=TRUE, an \code{os_iv_ci} confidence-set object.}
\seealso{\code{\link{iv.ci}}, \code{\link{iv.test}}}
\examples{
sens.ci.iv(rep(1,6),rep(0,6),2:7,rep(0,6),interval=c(-10,15),tol=0.1)
}
