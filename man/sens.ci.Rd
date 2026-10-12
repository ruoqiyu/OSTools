\name{sens.ci}
\alias{sens.ci}
\title{
Point estimate and confidence interval for sensitivity analysis in observational studies.
}
\usage{
sens.ci(y,Gamma=1,method='m',inner=0,trim=3,lambda=1/2,
    weight.par=c(1,1,1),alpha=0.05,alternative='two-sided',
    tol=NULL, interval=NULL)
}

\arguments{
\item{y}{
An n by J matrix, where (i) the rows are n matched sets, (ii) the first column is the treated response in a set, columns 2 to J contain the responses of controls in the same matched set.  Every set must have J-1 controls, and NAs are not allowed in y.
}
\item{Gamma}{
The sensitivity parameter: Gamma=1 for a randomization test, Gamma>1 for sensitivity bounds. This parameter gamma is denoted by the upper case Greek letter gamma in the literature, e.g., Rosenbaum (2007, 2014).
}
\item{method}{The method argument is the same as for the sens.analysis function.  See the help file for sens.analysis for more information.
}
\item{inner}{Inner trimming to increase design sensitivity.  See the help file for sens.analysis for more information.}
\item{trim}{Outer trimming for resistance to outliers.  See the help file for sens.analysis for more information.}
\item{lambda}{Observations are scaled by the lambda quantile of the absolute pair differences.  See the help file for sens.analysis for more information.}
\item{weight.par}{A vector of three parameters (m,m1,m2) that determine the weights, as in expression (5) in Rosenbaum (2014) or expression (8) in Rosenbaum (2011) where they are called m (for m), underline(m) for m1, and overline(m) for m2.
}
\item{alpha}{1-alpha is the coverage of the confidence interval.}
\item{alternative}{
If alternative='greater', the one-sided alternative hypothesis is that treatment increases the level of response by at least tau. If alternative='less', the one-sided alternative that the treatment decreases the level of response by at least tau. If alternative='two-sided', the alternative hypothesis is two-sided.
}
\item{tol}{For M-statistics and aligned ranks, the positive finite root-finding tolerance. For signed ranks, the positive finite grid spacing, which must not exceed the search width. If NULL, the interval width divided by 500000 is used.}
\item{interval}{Two finite, increasing search endpoints. If NULL, use the range from the smallest treated minus largest control response to the largest treated minus smallest control response. The required finite confidence bounds must lie within this range; widen it if a bound is not bracketed.}
}

\value{
\item{PointEstimate}{An interval of point estimates allowing for a bias of gamma in treatment assignment.}
\item{Confidence.Interval}{A confidence interval allowing for a bias of Gamma in treatment assignment. Signed-rank endpoints approximate the retained set on the supplied grid. One-sided intervals have an infinite endpoint on the untested side.}
}

\description{
Obtains point estimates and confidence intervals by inverting the hypothesis test in sens.analysis.  For instance, the confidence interval is the set of hypotheses not rejected by the test.  Estimates an additive treatment effect tau.  In a sensitivity analysis with gamma>1, one obtains an interval of point estimates consistent with a bias of at most gamma, and an interval containing all of the confidence intervals that are consistent with a bias of at most gamma.  
}

\examples{
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
X<-cbind(female,age,black,education,povertyr,bmi)
p<-glm(z~female+age+black+education+povertyr+bmi,family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
dist<-maha.dense(d$z,X)
dist<-add.caliper(dist, d$z, d$p, c(-.5,.15), stdev=TRUE, penalty=1000)
o<-match.os(d$z, dist, d, fine=d$education, ncontrol=1)
md<-o$data
sens.ci(cbind(md$homocysteine[md$z==1],md$homocysteine[md$z==0]),Gamma=1.2)
}


\references{
Huber, P. (1981) Robust Statistics.  New York: Wiley, 1981.

Maritz, J. S. (1979) Exact robust confidence intervals for location.  Biometrika 1979, 66, 163-166.

Rosenbaum, P. R. (2007) Sensitivity analysis for m-estimates, tests and confidence intervals in matched observational studies.  Biometrics, 2007, 63, 456-464. 

Rosenbaum, P. R. (2013) Impact of multiple matched controls on design sensitivity in observational studies.  Biometrics, 2013, 69, 118-127.

Rosenbaum, P. R. (2014) Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls.  Journal of the American Statistical Association, 109(507), 1145-1158

Rosenbaum, P. R. (2015). Two R packages for sensitivity analysis in observational studies. Observational Studies, 1(2), 1-17.  The Observational Studies journal is available free on-line.

}

\note{
Signed-rank calculations stop with an explanatory error if a required finite bound is truncated by the search interval, or if no effects are retained on the grid. Widen the interval or reduce the grid spacing as indicated. A finite bound may not be estimable for some sample sizes, significance levels or sensitivity parameters; increasing the range does not guarantee success. One-sided calculations check only their finite endpoint. Critical values use a direct upper-tail quantile to avoid subtracting tiny probabilities from one.

For M-statistics and aligned ranks, required roots must be bracketed by the finite search endpoints; otherwise a descriptive error identifies the missing bound. Nonfinite inputs, invalid significance levels and tuning arguments are rejected before searching. Redescending weights are not supported for confidence intervals. Existing test approximations remain unchanged.
}
