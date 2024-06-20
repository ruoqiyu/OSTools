\name{sens.ci}
\alias{sens.ci}
\title{
Point estimate and confidence interval for sensitivity analysis in observational studies.
}
\usage{
sens.ci(y,gamma=1,method='m',inner=0,trim=3,lambda=1/2,
    weight.par=c(1,1,1),alpha=0.05,alternative='two-sided',
    tol=NULL, interval=NULL)
}

\arguments{
\item{y}{
An n by J matrix, where (i) the rows are n matched sets, (ii) the first column is the treated response in a set, columns 2 to J contain the responses of controls in the same matched set.  Every set must have J-1 controls, and NAs are not allowed in y.
}
\item{gamma}{
The sensitivity parameter: gamma=1 for a randomization test, gamma>1 for sensitivity bounds. This parameter gamma is denoted by the upper case Greek letter gamma in the literature, e.g., Rosenbaum (2007, 2014).
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
\item{tol}{The sens.ci function calls the R function uniroot, and tol is the tol (or tolerance) parameter in that call. If tol=NULL, sens.ci picks a reasonable tolerance.}
\item{interval}{The sens.ci function calls the R function uniroot, and interval is the interval parameter in that call. If interval=NULL, sens.ci picks a reasonable interval.}
}

\value{
\item{PointEstimate}{An interval of point estimates allowing for a bias of gamma in treatment assignment.}
\item{CI}{An confidence interval allowing for a bias of gamma in treatment assignment.}
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
sens.ci(cbind(md$homocysteine[md$z==1],md$homocysteine[md$z==0]),gamma=1.2)
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
sens.ci inverts a test to obtain confidence intervals and point estimats; so, it calls senmw many times, solving several equations, and sens.ci is much slower than a single call to senmw. sens.ci finds point estimates and confidence intervals by searching for a value of the parameter tau in "interval" determining the solution tau.hat to an estimating equation with an error of "tol" in solving the equation.  If interval=NULL and tol=NULL, sens.ci tries to pick a reasonable finite interval and tol>0.  If concerned about these "reasonable values", set detail=TRUE, make the interval longer, the tol smaller, and wait longer for program to run.

Unlike sens.analysis, senms.ci does not permit redescending rank scores.
}

