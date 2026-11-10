\name{sens.analysis.tobit}
\alias{sens.analysis.tobit}
\alias{print.os_tobit_test}
\title{Matched-set sensitivity test for a censored reduction in a bounded response}
\usage{
sens.analysis.tobit(y, z, mset, Gamma = 1, tau = 0,
  alternative = c("two-sided", "greater", "less"),
  method = c("mean", "m"), inner = 0, trim = 3, lambda = 1/2)
\method{print}{os_tobit_test}(x, ...)
}
\arguments{
\item{x}{An object returned by the corresponding function.}
\item{...}{Additional print arguments (currently unused).}
\item{y}{Finite numeric work fractions in [0,1].}
\item{z}{Zero-one indicators, with exactly one treated individual per matched set.}
\item{mset}{Nonmissing matched-set labels; every set needs at least one control.}
\item{Gamma}{Finite sensitivity parameter at least one, for differential assignment.}
\item{tau}{Hypothesized reduction in [0,1]. Positive tau decreases the response.}
\item{alternative}{Direction for the reduction parameter, not for the raw treated-minus-control outcome: greater tests a larger reduction.}
\item{method}{Untrimmed mean scores or robust M-scores.}
\item{inner,trim,lambda}{M-score tuning: 0 <= inner <= trim, trim > 0, and 0 < lambda < 1. The mean method ignores these tuning values after validation.}
}
\description{Tests the model \eqn{Y(1)=\max(0,Y(0)-\tau)} using adjusted outcomes and a separable Normal approximation to matched-set sensitivity analysis. It is not a normal-error Tobit regression.}
\details{
For a candidate tau the adjusted response is pmax(0,y-(1-z)*tau).
For each possible treated individual in a set, its score is the sum of
within-set pair differences (or their odd M-score transforms), divided by
set size. These are the full-matching weights used by sensitivityfull;
with a common number of controls the mean statistic is proportional to the
usual treated-minus-control mean. Unequal-size sets do not receive equal
weights on their treated-minus-control means.

Mean scores use untrimmed differences. M-scores use the lambda quantile of
absolute ordered within-set differences (excluding self-differences). If
that quantile is zero but positive differences exist, the smallest positive
difference is used as scale and scale.fallback is TRUE. If all scores are
zero, the uninformative test returns p.value=1, including for a two-sided
test. The convention handles complete censoring and complete ties.

A greater reduction reverses the raw outcome direction. The two-sided
p-value doubles the smaller one-sided upper bound, capped at one.
P-values use the Normal upper tail directly, with a log-scale value also
returned. The reference is an asymptotic separable approximation, not an
exact finite-sample randomization distribution. The M-score fallback and
statistic weights must be considered when comparing published analyses.
}
\value{An object of class os_tobit_test with p.value (also pval), log.p.value,
deviate, statistic, expectation, variance, tau, Gamma, alternative, method,
adjusted outcomes in input order, scale, scale.fallback, direction-specific
centered scores, and call.}
\examples{
y <- rep(c(.3,.5,.5),12)
z <- rep(c(1,0,0),12)
mset <- rep(1:12,each=3)
sens.analysis.tobit(y,z,mset,tau=.2,alternative="greater")
}
\references{
Zubizarreta, J. R., Small, D. S. and Rosenbaum, P. R. (2014).
Isolation in the construction of natural experiments.
The Annals of Applied Statistics 8, 2096--2121. doi:10.1214/14-AOAS770.

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests and
confidence intervals in matched observational studies. Biometrics 63, 456--464.
}
\seealso{\code{\link{sens.ci.tobit}}, \code{\link{sens.analysis.full}}}
