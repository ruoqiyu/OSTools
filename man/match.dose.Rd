\name{match.dose}
\alias{match.dose}
\title{Nonbipartite Matching for Treatment Doses}
\description{Constructs disjoint pairs or full matched sets of subjects with similar covariates and separated treatment doses. Full matching uses a doubled-graph edge-cover construction and returns sets shaped as stars, including pairs.}
\usage{
match.dose(X, dose, method=c("full", "pair"),
  distance=c("rank", "mahalanobis", "euclidean"), lambda=0,
  near.exact=NULL, eta=1, mindist=0.01, penalty=100000, ndiscard=0)
}
\arguments{
\item{X}{Numeric covariate matrix, numeric data frame, or a vector representing one covariate. At least two subjects and one covariate are required. Values must be finite; encode categorical variables numerically and handle missing values before calling this function. Rows remain in their original order.}
\item{dose}{A finite numeric vector of treatment doses, one per row of X. For an ordinal exposure, supply numeric scores with the desired spacing.}
\item{method}{\code{"full"} retains all subjects in sets of at least two; \code{"pair"} forms disjoint pairs and optionally discards subjects.}
\item{distance}{\code{"rank"} uses the robust squared rank-based Mahalanobis distance from \code{maha.all}; \code{"mahalanobis"} uses squared ordinary Mahalanobis distance with a generalized inverse; \code{"euclidean"} uses ordinary Euclidean distance. Constant covariates contribute zero.}
\item{lambda}{Nonnegative finite penalty for each subject beyond the first two in a full matched set. Larger values favor more pairs. Must be zero for pair matching. Measured in the same units as the selected distance.}
\item{near.exact}{Optional column names or column numbers in X. For each selected covariate, an edge joining unequal values receives an extra cost of eta times the largest unpenalized distance. This is a pairwise mismatch penalty, not an aggregate fine-balance constraint.}
\item{eta}{Nonnegative finite multiplier for each near-exact mismatch penalty.}
\item{mindist}{Nonnegative minimum preferred absolute dose difference along each selected connection. Doses closer than this receive penalty. Equality with mindist is allowed.}
\item{penalty}{Nonnegative additive cost for a connection with dose difference below mindist. A finite value is a soft restriction and a warning reports selected violations. Use Inf to forbid such connections; an error occurs if the solver cannot return an admissible match.}
\item{ndiscard}{Number of subjects to discard in pair matching, an integer between zero and nrow(X)-2. The retained count must be even. Full matching requires zero. For an odd sample size, pair matching requires an odd ndiscard.}
}
\details{
For full matching, the original graph is duplicated, with cost \eqn{2(m_i+\lambda)} between subject i and its copy, where \eqn{m_i} is the cheapest admissible connection incident to i. Other cross-copy connections are forbidden. Minimum-cost nonbipartite pair matching on this graph is converted to an edge cover. Redundant connections are removed so each component is a star. A star contains a center joined to every other subject; no binary treated/control composition is implied.

The full-match objective is the sum of selected edge costs plus lambda times the sum of (set size minus two). The pair objective is the sum of selected edge costs. Pairwise costs include dose and near-exact penalties. The solver discretizes costs internally with \code{precision=8}; numerically tied optima may produce different equally good matches. Inadmissible returned connections are rejected.

All distance and penalty parameters use the native units of the chosen metric. Unlike the course helper file, this function does not multiply distances by 100 before constructing the graph. Existing tuning values from that file should be reconsidered.

The dose restriction concerns the selected star connections, not every pair of subjects within a full matched set. Different leaves may have identical doses. Distances are stored densely, with quadratic memory use; the doubled full-match graph has 2n vertices.
}
\value{A list with:
\item{mset}{Integer matched-set IDs, one per original input row. Discarded subjects have NA.}
\item{dose}{Input doses, in original order.}
\item{edges}{Data frame with original row positions i and j, set ID mset, unpenalized distance, penalized cost, and absolute dose.diff for each selected connection.}
\item{diagnostics}{One row per matched set, containing mset, size, center (original row position), mean.dist (mean over all distinct within-set pairs), star.dist (mean over selected connections), dose.diff (mean absolute dose difference over those connections), min.dose.diff, n.close (number below mindist), and the sum of edge costs. For pairs, the smaller row position is reported as the center.}
\item{mean.distance.before}{Mean unpenalized distance over all distinct input pairs.}
\item{discarded}{Original row positions of discarded subjects.}
\item{method, distance}{The matching method and distance metric used.}
\item{objective}{Selected matching objective, including lambda for full matching.}
\item{call}{The function call.}
}
\seealso{\code{\link{balance.check.dose}}, \code{\link{maha.all}}, \code{\link[nbpMatching]{nonbimatch}}}
\examples{
X <- cbind(age=c(20,21,30,31,45,46), income=c(3,4,5,6,7,8))
dose <- c(0,2,1,4,2,6)
fit <- match.dose(X, dose, method="full", lambda=1)
fit$mset
fit$diagnostics
balance.check.dose(X, dose, fit$mset)
pairs <- match.dose(X, dose, method="pair")
pairs$edges
}
