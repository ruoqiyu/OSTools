\name{sc.placebo}
\alias{sc.placebo}
\alias{print.os_sc_placebo}
\title{Descriptive In-Space Synthetic-Control Placebos}
\description{Applies the same outcome-only synthetic-control fitting rule to every candidate unit in turn, using all remaining units as donors. Reports post/pre RMSPE ratios and the treated unit's upper-tail rank fraction.}
\usage{
sc.placebo(outcomes, treated, pre, post=NULL, ridge=1e-8)
\method{print}{os_sc_placebo}(x, ...)
}
\arguments{
\item{outcomes}{Finite numeric matrix or all-numeric data frame, with chronologically ordered times in rows and candidate units in columns, including the actual treated unit. At least two rows and two columns are required. Rows must already align across units. Column names must be nonempty and unique if supplied; otherwise names \code{unit1}, \code{unit2}, etc. are generated.}
\item{treated}{One candidate column name or positive integer column index.}
\item{pre, post}{Fitting and evaluation row positions, specified as logical vectors or unique positive integer row indices. The default \code{post} selects all rows outside \code{pre}. Both selections must be nonempty and disjoint, and every pre row must precede every post row. See \code{\link{sc.fit}}.}
\item{ridge}{Numerical ridge used unchanged in every \code{sc.fit} call; default \code{1e-8}.}
\item{x}{An object returned by \code{sc.placebo}.}
\item{...}{Additional arguments passed to printing the rank summary.}
}
\details{
For each candidate column, fit using all other columns and compute the ratio of post RMSPE to pre RMSPE. The actual treated unit is included as a donor when fitting another candidate, just as in Chapter 17's symmetric relabeling exercise. No candidates are silently excluded for poor fit. A positive post RMSPE divided by zero pre RMSPE is \code{Inf}; a zero/zero ratio is undefined. If any candidate has an undefined ratio, \code{rank} and \code{rank.fraction} are \code{NA}, with \code{rank.defined=FALSE}. The full reference size is retained.

The rank counts all candidates, including the treated unit, whose ratio is at least as large as the treated unit's. Ties are included. The rank fraction is that count divided by the number of candidates. This function deliberately supplies no \code{p.value} field: the fraction is descriptive by default. Under a sharp null, a uniformly assigned treated label over these candidates and an identical, symmetrically applied analysis rule give the usual randomization interpretation. That assignment model must be justified separately; an observed policy adoption is not automatically a lottery. Unequal assignment probabilities, asymmetric donor restrictions, outcome-dependent candidate selection or fitted tuning choices can invalidate that interpretation. Zero/zero ratios would need a prespecified symmetric statistic convention before a randomization test could be defined.

Weights use only \code{pre} outcomes. Ratios with a nearly zero denominator can be very sensitive to floating-point error or small data changes; inspect the separate RMSPE values rather than interpreting a large ratio alone. There is no automatic inference for the average treatment effect, and a placebo rank does not establish the validity of counterfactual prediction. The fit and donor assumptions of \code{sc.fit} still apply. In-time placebos and donor-omission checks can be computed by calling \code{sc.fit} with appropriate rows or columns instead.
}
\value{An object of class \code{os_sc_placebo}. Fields are \code{treated} (unit label), \code{statistics} (a data frame with unit, pre/post RMSPE, ratio and post mean gap), named \code{ratios}, \code{rank}, \code{rank.fraction}, \code{reference.size}, \code{rank.defined}, \code{weights} (donors in rows, candidate treated units in columns; diagonal zero), \code{gaps} (times in rows, candidates in columns), \code{pre}, \code{post}, \code{ridge}, \code{interpretation} and \code{call}.}
\references{
Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program. \emph{Journal of the American Statistical Association}, 105, 493--505. doi:10.1198/jasa.2009.ap08746.
}
\seealso{\code{\link{sc.fit}}}
\examples{
outcomes <- cbind(Treated=c(2,3,4,5,9,10),
                  A=c(1,2,3,4,5,6), B=c(3,4,5,6,7,8))
placebo <- sc.placebo(outcomes,treated="Treated",pre=1:4)
placebo
placebo$statistics
}
