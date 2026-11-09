\name{iv.test}
\alias{iv.test}
\alias{print.os_iv_test}
\title{Rank Tests of Instrumental-Variable Adjusted Outcomes}
\description{Tests a candidate additive effect using adjusted outcomes, with rank-sum inference for two independent instrument groups or signed-rank inference for matched pairs. Paired analyses can also bound sensitivity to hidden bias in instrument assignment.}
\usage{
iv.test(y, treatment=NULL, instrument=NULL, mset=NULL, beta=0,
  method=c("signedrank", "ranksum"), Gamma=1,
  alternative=c("two.sided", "greater", "less"),
  exact=NULL, correct=TRUE, zero.method=c("wilcox", "pratt"),
  d=NULL, z=NULL)
\method{print}{os_iv_test}(x, ...)
}
\arguments{
\item{y, treatment}{Finite numeric or logical outcome and treatment-received vectors, denoted Y and Z in Chapter 15. For paired analysis, they can instead be numeric matrices of identical dimensions with two columns, the higher-instrument member first.}
\item{instrument}{Numeric or logical instrument vector, denoted A in Chapter 15. For rank-sum inference it must have exactly two levels; the larger level denotes encouragement. With matched-pair labels it can be continuous, but must differ within every retained pair. Omit for paired matrix inputs.}
\item{mset}{Matched-pair labels for vector inputs with \code{method="signedrank"}. Each retained label must occur twice. Missing labels identify unmatched subjects, who are excluded. Outcome, treatment and instrument inputs must still be finite for all supplied subjects. Omit for matrix inputs or rank-sum analysis.}
\item{d, z}{Compatibility aliases from version 0.1.10: d means treatment and z means instrument. Their old meanings are preserved. Supply only one name for each role; specifying a canonical argument and its alias together is an error.}
\item{beta}{One finite candidate effect.}
\item{method}{Signed-rank or rank-sum test.}
\item{Gamma}{Finite sensitivity parameter at least 1. Values above 1 require paired signed-rank analysis and use a normal approximation.}
\item{alternative}{Direction of the adjusted-outcome contrast for the higher-instrument group. This is not automatically the direction of an effect alternative if the first stage is negative.}
\item{exact}{Logical, or \code{NULL}. At Gamma 1, automatic exact inference uses untied samples below 50 observations per group, or below 50 nonzero differences for pairs without zeros. Ties, zeros and Pratt scores use a normal approximation; explicitly requesting exact inference then issues a warning. Exact sensitivity bounds for Gamma above 1 are not implemented.}
\item{correct}{Apply a 0.5 continuity correction in normal approximations. Ignored by exact tests.}
\item{zero.method}{For signed ranks, \code{"wilcox"} removes zero adjusted differences before ranking; \code{"pratt"} includes them when ranking absolute differences, but gives them no signed contribution. Pratt inference is a normal approximation.}
\item{x}{A returned test object.}
\item{...}{Additional printing arguments.}
}
\details{
Chapter 15 uses A for the instrument and Z for treatment received. The adjusted responses are \eqn{Y-\beta Z}. Use \code{treatment=Z, instrument=A}; do not pass the treatment to the legacy \code{z} keyword. Positional second and third arguments retain their original treatment/instrument order. Independent groups are compared using \code{stats::wilcox.test}; at Gamma 1 the default zero convention also uses that function for pairs. Paired sensitivity bounds use conditional signed-rank normal moments. \code{p.value} is the upper, worst-case sensitivity bound when Gamma exceeds 1. Zero adjusted differences throughout, or identical outcomes in the two independent groups, give p-value 1.

A causal interpretation requires appropriate instrument assignment, exclusion and an additive effect model. For binary treatment, a common effect among compliers suffices with the additional assumptions described in the textbook. These rank procedures do not generally identify the mean complier effect under unrestricted heterogeneity. In matched pairs, assignment of the higher instrument value must be as if random conditional on the pair; Gamma relaxes that assignment assumption, not exclusion. The functions do not test whether these assumptions hold.

Higher instrument values need not produce higher observed treatment values in every pair. Confidence sets should be obtained by inversion rather than by assuming a particular first-stage sign. Normal approximations with few informative pairs can be inaccurate.
}
\value{An object of classes \code{os_iv_test} and \code{htest}, containing \code{statistic}, \code{p.value}, \code{p.value.bounds} (lower and upper), \code{null.value}, \code{alternative}, \code{method}, \code{Gamma}, the observed \code{first.stage} difference, group sizes \code{n}, the count \code{n.discarded}, the actual \code{exact} flag, \code{correct}, \code{zero.method} and \code{call}.}
\seealso{\code{\link{iv.ci}}, \code{\link{sens.ci.iv}}, \code{\link[stats]{wilcox.test}}}
\examples{
y <- cbind(c(2,3,4,5,6,7), c(0,0,0,0,0,0))
Z <- cbind(rep(1,6), rep(0,6))
iv.test(y,treatment=Z,beta=0)
iv.test(y,treatment=Z,beta=0,Gamma=2,exact=FALSE,correct=FALSE)
iv.test(y=c(2,4,6,1,3,5), treatment=c(1,1,1,0,0,0),
        instrument=c(1,1,1,0,0,0), method="ranksum")
}
