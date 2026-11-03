\name{combine.pvalues}
\alias{combine.pvalues}
\title{Fisher and Truncated-Product Combinations of P-Values}
\usage{
combine.pvalues(p, method = c("fisher", "truncated"),
  trunc = 0.2, log.p = FALSE)
}
\arguments{
\item{p}{Nonempty numeric vector of component p-values in \eqn{[0,1]}, or their natural logarithms when \code{log.p = TRUE}. Missing values are not permitted. If supplied, names must be unique and nonempty.}
\item{method}{\code{"fisher"} for Fisher's combination or \code{"truncated"} for the truncated product.}
\item{trunc}{Prespecified cutoff in \eqn{(0,1]}. Fisher's method always uses cutoff 1.}
\item{log.p}{Whether the inputs are log probabilities. Log inputs may equal negative infinity, representing zero.}
}
\description{Combines independent p-values, or jointly valid evidence-factor bounds, to test the intersection null that all component null hypotheses are true.}
\details{
Fisher's statistic is \eqn{-2\sum_i\log p_i}, calibrated by a chi-squared distribution with twice the number of components as its degrees of freedom. The truncated product uses
\deqn{W_c=\prod_{i:p_i\leq c}p_i,}
with the empty product equal to 1. The cutoff must be fixed before examining the p-values. Setting \code{trunc = 1} gives Fisher's method.

Under independent uniform null p-values, the number of retained components has a binomial distribution. Conditional on this number, the negative log product after scaling by the cutoff has a gamma distribution. The implementation sums these nonnegative mixture terms on the log scale, without forming the product of small probabilities or subtracting a nearly unit lower-tail probability. It includes the point mass at product 1. This is the same calibration as the truncated-product formula, not Fisher's method applied after dropping large p-values.

\code{log.p.value} can remain finite when the ordinary p-value underflows to zero. Supplying logarithms avoids losing small input probabilities before the function is called. A zero supplied on the ordinary scale is treated as zero; its unrounded value cannot be recovered.

Validity requires independent null p-values or the joint stochastic bounds required for the monotone combination, as with appropriate evidence factors. Marginally valid p-values alone do not suffice under arbitrary dependence, and a small empirical correlation does not establish validity. The function cannot verify these assumptions. A rejection supports at least one false component null; it does not reject every component null. Sensitivity bounds may be supplied only when their joint validity is justified.
}
\value{An \code{htest} object with \code{statistic}, \code{parameter}, \code{p.value}, \code{log.p.value}, \code{alternative}, \code{method}, and \code{data.name}. It also contains \code{p.values} and \code{log.p.values} in input order, \code{combination}, the effective \code{trunc}, \code{n.retained}, and an \code{assumptions} statement. The truncated-product statistic is reported as its logarithm; the ordinary product need not be representable.}
\references{
Zaykin, D. V., Zhivotovsky, L. A., Westfall, P. H. and Weir, B. S. (2002). Truncated product method for combining P-values. \emph{Genetic Epidemiology}, 22, 170--185. \doi{10.1002/gepi.0042}.

Hsu, J. Y., Small, D. S. and Rosenbaum, P. R. (2013). Effect modification and design sensitivity in observational studies. \emph{Journal of the American Statistical Association}, 108, 135--148. \doi{10.1080/01621459.2012.742018}.

Rosenbaum, P. R. (2011). Some approximate evidence factors in observational studies. \emph{Journal of the American Statistical Association}, 106, 285--295. \doi{10.1198/jasa.2011.tm10422}.
}
\examples{
# Chapter 13: chromium evidence-factor bounds (rounded as in the text).
p <- c(workers = 0.0588, exposure = 0.0819)
combine.pvalues(p, method = "fisher")$p.value
combine.pvalues(p, method = "truncated", trunc = 0.2)$p.value

# The log result remains usable when the ordinary probability underflows.
combine.pvalues(c(-1000, -1200), log.p = TRUE)$log.p.value
}
\seealso{\code{\link{partial.conjunction}}, \code{\link{sens.analysis}}}
