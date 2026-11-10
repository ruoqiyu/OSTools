\name{sc.fit}
\alias{sc.fit}
\alias{print.os_sc_fit}
\title{Outcome-Only Synthetic Control}
\description{Fits a nonnegative, sum-to-one weighted combination of donor outcome histories using only specified pretreatment rows. Returns weights, a synthetic outcome path, gaps and fit summaries.}
\usage{
sc.fit(y, donors, pre, post=NULL, ridge=1e-8)
\method{print}{os_sc_fit}(x, ...)
}
\arguments{
\item{y}{Finite numeric vector of the treated unit's outcomes, ordered chronologically.}
\item{donors}{Finite numeric matrix or all-numeric data frame with time in rows and donors in columns. Its rows must already align with \code{y}; names are not used to align or sort observations. There must be at least two rows and one donor. Column names, if supplied, must be nonempty and unique; otherwise unit labels are generated.}
\item{pre}{Logical vector of length \code{length(y)}, or unique positive integer row indices, selecting at least one fitting period. These are row positions, not calendar dates.}
\item{post}{Evaluation row positions in the same format as \code{pre}. By default all rows outside \code{pre}. Must select at least one period, be disjoint from \code{pre}, and occur strictly after every fitting period. An explicit subset can leave other rows unused in the summaries, for example to evaluate an in-time placebo before the actual intervention.}
\item{ridge}{Finite nonnegative numerical ridge, default \code{1e-8}. It is added to the diagonal of the pretreatment donor cross-product matrix. A zero ridge requires a positive-definite cross-product; otherwise optimization fails explicitly. The ridge is scale dependent and should be kept fixed across comparisons on the same outcome scale.}
\item{x}{An object returned by \code{sc.fit}.}
\item{...}{Additional arguments passed to printing the numerical summary.}
}
\details{
The weights minimize
\deqn{\sum_{t\in\mathrm{pre}}(y_t-\sum_j w_j X_{tj})^2+\mathrm{ridge}\sum_j w_j^2}
subject to nonnegative weights summing to one. \pkg{quadprog} solves the quadratic program; tiny numerical constraint violations are clipped and the weights renormalized. The default ridge and normalization reproduce Chapter 17's original outcome-only helper. This is not the predictor-weighted \pkg{Synth} specification. It includes no intercept, predictor balancing or covariate adjustment.

Weights use only fitting rows. Predictions and gaps are returned for every supplied row, including any rows omitted from the pre/post summaries. No missing outcomes are silently discarded. The gap is observed minus synthetic outcome; its post-period mean is an equally weighted average over the selected dates. RMSPE is the square root of the mean squared gap over the indicated dates.

A causal interpretation requires credible donor eligibility, no anticipation or spillovers, and a stable relationship between the untreated treated-unit outcome and the weighted donors. Good pretreatment fit does not verify those conditions. Neither RMSPE nor the reported mean gap is a standard error or a confidence interval. All units must remain untreated during fitting, and donors during evaluation. Users are responsible for those substantive conditions and for the chronological alignment of the rows.
}
\value{An object of class \code{os_sc_fit}, with \code{weights} (named donor weights), \code{synthetic}, \code{gap}, \code{outcome}, \code{pre} and \code{post} (selected row indices), \code{pre.rmspe}, \code{post.rmspe}, \code{post.mean.gap}, \code{ridge} and \code{call}.}
\references{
Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program. \emph{Journal of the American Statistical Association}, 105, 493--505. doi:10.1198/jasa.2009.ap08746.
}
\seealso{\code{\link{sc.placebo}}, \code{\link[quadprog]{solve.QP}}}
\examples{
donors <- cbind(A=c(1,2,3,4,5,6), B=c(3,1,4,2,6,3))
y <- as.vector(donors \%*\% c(.25,.75))+c(0,0,0,0,2,2)
fit <- sc.fit(y,donors,pre=1:4)
fit
fit$weights
fit$gap
# Hold out an untreated date, omitting the true post-treatment dates.
sc.fit(y,donors,pre=1:3,post=4)
}
