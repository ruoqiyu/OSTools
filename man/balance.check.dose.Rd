\name{balance.check.dose}
\alias{balance.check.dose}
\title{Covariate Balance Within Dose-Matched Sets}
\description{Compares covariate means for the higher- and lower-dose halves of matched sets and reports standardized mean differences before and after matching.}
\usage{balance.check.dose(X, dose, mset)}
\arguments{
\item{X}{Numeric covariate matrix, numeric data frame, or vector for a single covariate. Include the entire original sample, in the same row order as dose and mset. At least four subjects and one covariate are required; all values must be finite.}
\item{dose}{Finite numeric vector of treatment doses.}
\item{mset}{Matched-set labels, one per row of X, such as \code{fit$mset} from \code{match.dose}. Labels may be numeric, character, or factors and need not be consecutive. NA marks discarded subjects. Every retained set must have at least two members.}
}
\details{
Before matching, all subjects are sorted by dose. The first floor(n/2) form the lower-dose group and the remaining subjects form the higher-dose group. After matching, the same rule is applied separately in each retained matched set. An odd-sized set assigns its extra subject to the higher-dose group. Dose ties are broken deterministically by original row order. A warning identifies sets with no dose variation.

Within each matched set, covariate means are computed separately for its higher- and lower-dose groups. The reported after-matching means average these set-specific means with equal weight for each set, following the course helper's convention. They are not pooled subject-weighted means when set sizes differ.

Both before- and after-matching differences are divided by the same full-sample reference standard deviation: the square root of the average of the variances in the original higher- and lower-dose groups. Thus discarded subjects still contribute to the before-matching comparison and reference scale. Four original subjects are required so each reference group has at least two observations.

If the reference scale and mean difference are both zero, the SMD is zero. A nonzero difference with zero reference scale produces a signed infinite SMD. The sign is always higher dose minus lower dose. This diagnostic describes the chosen within-set split; it is not a test of a dose-response effect or a guarantee of balance across all dose levels.
}
\value{A numeric matrix with one row per covariate and columns \code{Before Higher Mean}, \code{Before Lower Mean}, \code{After Higher Mean}, \code{After Lower Mean}, \code{Before SMD}, and \code{After SMD}. Unnamed covariates are labeled X1, X2, and so on.}
\seealso{\code{\link{match.dose}}, \code{\link{love.plot}}, \code{\link{balance.check}}}
\examples{
X <- cbind(age=c(20,21,30,31,45,46), income=c(3,4,5,6,7,8))
dose <- c(0,2,1,4,2,6)
fit <- match.dose(X, dose, method="full", lambda=1)
bal <- balance.check.dose(X, dose, fit$mset)
bal
love.plot(abs(bal[,"Before SMD"]), abs(bal[,"After SMD"]), rownames(bal))
}
