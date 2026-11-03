\name{closed.testing}
\alias{closed.testing}
\alias{print.os_closed_testing}
\title{Closed Testing of Component No-Effect Hypotheses}
\usage{
closed.testing(p, method = c("truncated", "fisher"), trunc = 0.2,
  alpha = 0.05, log.p = FALSE, max.groups = 15L)
\method{print}{os_closed_testing}(x, digits = 4, ...)
}
\arguments{
\item{p}{Nonempty numeric vector of component p-values or upper bounds. Names, if supplied, must be unique and nonempty.}
\item{method}{Local combination used for intersections of two or more hypotheses.}
\item{trunc}{Prespecified truncated-product cutoff in \eqn{(0,1]}. Fisher uses cutoff 1.}
\item{alpha}{Significance level strictly between zero and one. Rejection uses p-value at most \code{alpha}.}
\item{log.p}{Whether \code{p} contains natural logarithms. Negative infinity represents zero.}
\item{max.groups}{Maximum permitted number of component hypotheses, an integer from 1 to 20, default 15. Closure requires \eqn{2^G-1} local tests.}
\item{x}{Object returned by \code{closed.testing}.}
\item{digits}{Printing precision.}
\item{...}{Additional arguments passed to data-frame printing.}
}
\description{Enumerates the nonempty intersections of component hypotheses and applies the closure principle. Returns adjusted p-values for individual components and an explicit table of the local intersection tests.}
\details{
For a subset \eqn{K}, a local test combines its component p-values using \code{\link{combine.pvalues}}. A singleton always uses its original p-value, even if it exceeds the truncation cutoff. For component \eqn{j}, the adjusted p-value is the maximum local p-value over all subsets containing \eqn{j}. The global-null test is only one of those tests; its rejection need not identify an individual component. The intersection table's \code{local.reject} denotes a local rejection, not rejection by closure.

Closure controls the familywise error rate when every relevant intersection has a valid local level-alpha test. Fisher and truncated-product calibration require independent true-null p-values or the necessary joint stochastic bounds; marginal validity alone under arbitrary dependence is insufficient. The method and cutoff must be fixed independently of the signed results. With data-dependent groups, the entire discovery procedure must satisfy the corresponding conditional validity argument. These conditions cannot be checked from \code{p}.

Combinations, maxima and decisions are computed on the log scale. Extremely small probabilities can underflow in the ordinary-scale output while their logarithms remain available. The guard on group count prevents accidental exponential enumeration; increasing it is an explicit computational choice.
}
\value{An \code{os_closed_testing} list. \code{tests} contains original component labels, \code{p.value}, \code{p.adjusted}, \code{reject}, \code{log.p.value}, and \code{log.p.adjusted}. \code{intersections} contains subset IDs, display labels, sizes, local p-values and log p-values, \code{local.reject}, and an \code{indices} list column identifying each subset without relying on display labels. Other fields include \code{global.p.value}, \code{global.log.p.value}, \code{alpha}, \code{combination}, \code{trunc} and \code{assumptions}.}
\references{
Marcus, R., Peritz, E. and Gabriel, K. R. (1976). On closed testing procedures with special reference to ordered analysis of variance. \emph{Biometrika}, 63, 655--660. \doi{10.1093/biomet/63.3.655}.

Hsu, J. Y., Zubizarreta, J. R., Small, D. S. and Rosenbaum, P. R. (2015). Strong control of the familywise error rate in observational studies that discover effect modification by exploratory methods. \emph{Biometrika}, 102, 767--782. \doi{10.1093/biomet/asv034}.
}
\examples{
closed.testing(c(G1 = .01, G2 = .30, G3 = .04))
closed.testing(c(A = -1000, B = log(.3)), log.p = TRUE)$tests
}
\seealso{\code{\link{combine.pvalues}}, \code{\link{partial.conjunction}}, \code{\link{sens.analysis.groups}}}
