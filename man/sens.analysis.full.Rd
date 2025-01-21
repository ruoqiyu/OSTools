\name{sens.ci.full}
\alias{sens.ci.full}
\title{
Sensitivity analysis in full matching.
}
\usage{
sens.anallysis.full(y,z,mset,Gamma=1,method='m',alternative='two-sided',
tau=0,inner=0,trim=3,lambda=1/2)
}

\arguments{
\item{y}{A vector of outcomes in the matched data.}
\item{z}{A vector of treatment indicators in the matched data.}
\item{mset}{A vector of matched set information.}
\item{Gamma}{
The sensitivity parameter: Gamma=1 for a randomization test, Gamma>1 for sensitivity bounds. This parameter Gamma is denoted by the upper case Greek letter Gamma in the literature, e.g., Rosenbaum (2007, 2014).
}
\item{method}{Choose method='mh' for a binary outcome and method='m' for a continuous outcome.}
\item{alternative}{
If alternative='greater', the one-sided alternative hypothesis is that treatment increases the level of response by at least tau. If alternative='less', the one-sided alternative that the treatment decreases the level of response by at least tau. If alternative='two-sided', the alternative hypothesis is two-sided.
}
\item{tau}{
If tau=0, test the null hypothesis of no treatment effect.  If tau is not 0, test the null hypothesis that the treatment effect is an additive shift of tau against the alternative that the effect is larger than tau in the direction specified by alternative.
}
\item{inner}{Inner trimming to increase design sensitivity.  See the help file for sens.analysis for more information.}
\item{trim}{Outer trimming for resistance to outliers.  See the help file for sens.analysis for more information.}
\item{lambda}{Observations are scaled by the lambda quantile of the absolute pair differences.  See the help file for sens.analysis for more information.}
}

\value{
\item{pval}{The upper bound on the 1-sided P-value.}
\item{data.array}{A 3-dimensional contingency table in array form where each dimension is at least 2 for Mantel-Haenszel test.}
\item{deviate}{The deviate that was compared to the Normal distribution to produce pval.}
\item{statistic}{The value of the M-statistic.}
\item{expectation}{The maximum expectation of the M-statistic for the given Gamma.}
\item{variance}{The maximum variance of the M-statistic among treatment assignments that achieve the maximum expectation. Part of the separable approximation.}
}

\description{
Obtains point estimates and confidence intervals for m test in full matching.
}

\examples{
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
X<-cbind(female,age,black,education,povertyr,bmi)
p<-glm(z~female+age+black+education+povertyr+bmi,family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
dist<-maha.dense(d$z,X,matrix=TRUE)
dist<-add.caliper(dist,d$z,d$p,c(-.2,.2),stdev=TRUE,penalty=1000)
o<-match.os.full(d$z,dist,d)
md<-o$data
sens.analysis.full(md$homocysteine,md$z,md$mset,Gamma=1.2,method='m')
}


\references{
Huber, P. (1981) Robust Statistics.  New York: Wiley, 1981.

Maritz, J. S. (1979) Exact robust confidence intervals for location.  Biometrika 1979, 66, 163-166.

Rosenbaum, P. R. (2007) Sensitivity analysis for m-estimates, tests and confidence intervals in matched observational studies.  Biometrics, 2007, 63, 456-464. 

Rosenbaum, P. R. (2013) Impact of multiple matched controls on design sensitivity in observational studies.  Biometrics, 2013, 69, 118-127.

Rosenbaum, P. R. (2014) Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls.  Journal of the American Statistical Association, 109(507), 1145-1158

}


