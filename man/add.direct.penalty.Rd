\name{add.direct.penalty}
\alias{add.direct.penalty}
\title{Add a directional penalty to a distance object}
\usage{
add.direct.penalty(dist, z, dx, positive = TRUE, penalty = 1)
}
\arguments{
\item{dist}{  A distance object with three components: d, start, end, typically created by maha_dense or maha_sparse. d[i] gives the distance between the (start[i])th treated and the (end[i]-sum(z))th control.}
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.  Must have
treated subjects (z=1) before controls (z=0).}
\item{dx}{ A vector of with length(z)=length(dx) giving the variable used to define the caliper.  For instance,
dx might be the propensity score.}
\item{positive}{ If positive = TRUE, a treated-minus-control difference in dx that is positive is
increased by penalty, but if positive = FALSE a a treated-minus-control difference in dx that is negative is
increased by penalty.  Zero differences are never penalized.}
\item{penalty}{ The number added to a distance when the desired direction is violated. }
}
\value{ Returns a new distance matrix that is the sum of dmat and the penalties for direction violations.}
\description{
Add a directional penalty to a distance object.

}
\examples{
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
X<-cbind(female,age,black,education,povertyr,bmi)
p<-glm(z~female+age+black+education+povertyr+bmi,family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
dist<-maha.dense(d$z,X)
head(dist$d)
dist<-add.direct.penalty(dist, d$z, d$p, positive=TRUE, penalty = 1)
head(dist$d)
}
\references{

Yu, R., & Rosenbaum, P. R. (2019). Directional penalties for optimal matching in observationalstudies.Biometrics, 75(4), 1380-1390.

}
