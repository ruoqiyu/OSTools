\name{dilated.effect.matchedpair}
\alias{dilated.effect.matchedpair}
\title{Randomization inference for the dilated treatment effect for matched pairs}
\usage{
dilated.effect.matchedpair(treated.r,control.r,k,Delta0=0,alternative="greater",
alpha=0.05,range=c(-100,100),grid=seq(0,10000,1),tol=0.0001)
}
\arguments{
\item{treated.r}{ A vector of responses in the treated group.}
\item{control.r}{ A vector of responses in the control group.}
\item{k}{ Quantile of interest for the potential outcome under control.}
\item{Delta0}{ Dilated treatment effect value in the null hypothesis.}
\item{alternative}{ Direction of the alternative hypothesis. It can be one of 'greater' and 'less'.}
\item{alpha}{ Coverage probability of confidence interval is 1-alpha.}
\item{range}{ Range to search for the root in the uniroot function, which is used to determine the two end points of the confidence interval.}
\item{grid}{ Consider grid of values, find smallest value such that test statistic is less than its expectation, and largest value such that test statistic is greater than its expectation, and average these values.}
\item{tol}{ Tolerance to estimate the Hodges-Lehmann point estimate.}
}
\value{
\item{pvalue}{P-value for testing Delta(rho)=Delta0 against the specified alternative hypothesis.}
\item{conf.interval}{Lower bound and upper bound of the (1-alpha) confidence interval.}
\item{hl.estimate}{The Hodges-Lehmann point estimate.}
}
\description{
Construct hypothesis testing of dilated treatment effect for matched pairs, search for endpoints of lower and upper bounds for confidence intervals, and estimate the Hodges-Lehmann point estimate.
}
\examples{
cadmium=c(107143,33679,18836,173,389,1144,513,211,24288,67632,488,700,
          328,98,122,2302,10208,892,2803,201,148,522,941)
hospital=c(311,338,159,110,226,305,222,242,250,256,135,96,
           142,120,376,173,178,213,257,81,199,114,247)

k=23
dilated.effect.matchedpair(cadmium,hospital,k,Delta0=100,
                           range=c(-10^3,10^3),grid=seq(0,10^3,1))
}

