\name{sens.ci.full}
\alias{sens.ci.full}
\title{
Point estimate and confidence interval for sensitivity analysis in full matching.
}
\usage{
sens.ci.full(y,z,mset,Gamma=1,inner=0,trim=3,lambda=1/2,alpha=0.05,
alternative='two-sided')
}

\arguments{
\item{y}{A vector of outcomes in the matched data.}
\item{z}{A vector of treatment indicators in the matched data.}
\item{mset}{A vector of matched set information.}
\item{Gamma}{
The sensitivity parameter: Gamma=1 for a randomization test, Gamma>1 for sensitivity bounds. This parameter Gamma is denoted by the upper case Greek letter Gamma in the literature, e.g., Rosenbaum (2007, 2014).
}
\item{inner}{Inner trimming to increase design sensitivity.  See the help file for sens.analysis for more information.}
\item{trim}{Outer trimming for resistance to outliers.  See the help file for sens.analysis for more information.}
\item{lambda}{Observations are scaled by the lambda quantile of the absolute pair differences.  See the help file for sens.analysis for more information.}

\item{alpha}{1-alpha is the coverage of the confidence interval.}
\item{alternative}{
If alternative='greater', the one-sided alternative hypothesis is that treatment increases the level of response by at least tau. If alternative='less', the one-sided alternative that the treatment decreases the level of response by at least tau. If alternative='two-sided', the alternative hypothesis is two-sided.
}
}

\value{
\item{PointEstimate}{An interval of point estimates allowing for a bias of Gamma in treatment assignment.}
\item{CI}{An confidence interval allowing for a bias of Gamma in treatment assignment.}
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
sens.ci.full(md$homocysteine,md$z,md$mset,Gamma=1.2)
}


\references{
Huber, P. (1981) Robust Statistics.  New York: Wiley, 1981.

Maritz, J. S. (1979) Exact robust confidence intervals for location.  Biometrika 1979, 66, 163-166.

Rosenbaum, P. R. (2007) Sensitivity analysis for m-estimates, tests and confidence intervals in matched observational studies.  Biometrics, 2007, 63, 456-464. 

Rosenbaum, P. R. (2013) Impact of multiple matched controls on design sensitivity in observational studies.  Biometrics, 2013, 69, 118-127.

Rosenbaum, P. R. (2014) Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls.  Journal of the American Statistical Association, 109(507), 1145-1158

}


