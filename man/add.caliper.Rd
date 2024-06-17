\name{add.caliper}
\alias{add.caliper}
\title{Add a caliper, that need not be symmetric, to a distance object.}
\usage{
add.caliper(dist, z, dx, rg, stdev = FALSE, penalty = 1000, constant = TRUE)
}
\arguments{
\item{dist}{ A distance matrix or a distance object with three components: d, start, end, typically created by maha_dense or maha_sparse. d[i] gives the distance between the (start[i])th treated and the (end[i]-sum(z))th control.}
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.  Must have
treated subjects (z=1) before controls (z=0).}
\item{dx}{ A vector of with length(z)=length(dx) giving the variable used to define the caliper.  For instance,
dx might be the propensity score.}
\item{rg}{ A vector with length(rg)=2 such that rg[1] <= 0 <= rg[2].  If the treated-minus-control
difference in dx is < rg[1] or > rg[2], then penalty is added to the distance.  If treated individuals
have dx higher than controls, then you want to set rg[2] < -rg[1], so that you tolerate smaller
positive differences and larger negative differences. If the user only specify one value for rg, assume the range is symmetric that rg[1]=-abs(rg) and rg[2]=abs(rg).}
\item{stdev}{ If stdev = TRUE, rg is interpreted in units of an equally weighted pooled
standard deviation; that is,rg is replaced by rg*sp where sp is
sqrt((var(dx[z==1])+var(dx[z==0]))/2).}
\item{penalty}{ The number added to a distance when the caliper is violated.  A large penalty,
like the default value of penalty = 1000, will try to enforce the caliper to the extent that this
is feasible.  Small penalties can slightly tilt a match in a desired direction. }
\item{constant}{ An indicator to denote whether the same penalty is used for all violations of the caliper or the penalty is proportional to the scale of deviations from the caliper. }
}
\value{ Returns a new distance object whose distance component d is updated by the sum of d and the penalties for caliper violations.}
\description{
Imposes a caliper, that need not be symmetric, on p using a penalty function, adding the penalty to a distance matrix dmat and returning a new distance matrix.

This symmetric version of this function is discussed in Rosenbaum (2010).

}
\examples{
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
X<-cbind(female,age,black,education,povertyr,bmi)
p<-glm(z~female+age+black+education+povertyr+bmi,family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
dist0<-maha.dense(d$z,X)
#symmetric caliper
dist1<-add.caliper(dist0, d$z, d$p, c(-.3,.3), stdev = TRUE, penalty = 1000)
head(dist1$d)
#asymmetric caliper
dist2<-add.caliper(dist0, d$z, d$p, c(-.5,.1), stdev = TRUE, penalty = 1000)
head(dist2$d)
}
\references{

Yu, R., & Rosenbaum, P. R. (2019). Directional penalties for optimal matching in observationalstudies.Biometrics, 75(4), 1380-1390.

}
