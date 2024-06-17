\name{love.plot}
\alias{love.plot}
\title{Graphical way of showing how the covariate balance measure changed from before matching to after.}
\usage{
love.plot(before,after,covariates,name='abs.stand.diff',
sizev=2,shapev=c(4,1),linex=c(0.1,0.2))
}
\arguments{
\item{before}{ A vector with length(covariates) values for the covariate balance measure for each covariate before matching.}
\item{after}{ A vector with length(covariates) values for the covariate balance measure for each covariate after matching.}
\item{covariates}{ A vector containing names of the covariates.}
\item{name}{ A charactor for the name of the covariate balance measure. The name is used as the x-axis label in the output plot.}
\item{sizev}{ A numerical value for the size of the points in the plot. The default value is 2.}
\item{shapev}{ A vector of two numerical values for the shapes of the points before and after matching.}
\item{linex}{ A vector of numerical values for the x intercept of vertical dashed lines. The default is (0.1,0.2) for the ideal cutoff and acceptable cutoff for absolute standardized mean differences.}
}
\value{ A figure produced by ggplot2.
}
\description{
The function is used to create a graphical illustration of how the covariate balance measure changed from before matching to after. Mostly commonly, this kind of plot (invented by Thomas Love) is used to compare the standardized mean differences before and after matching.
}
\examples{
\donttest{
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
X<-cbind(female,age,black,education,povertyr,bmi)
p<-glm(z~female+age+black+education+povertyr+bmi,family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
dist<-maha.dense(d$z,X)
o<-match.os(d$z, dist, d)
matcheddata<-o$data
Xm<-subset(matcheddata, select=c('female','age','black','education','povertyr','bmi','p'))
btb=balance.check(cbind(X,p),Xm,d$z,matcheddata$z)
love.plot(btb[,5],btb[,4],c('female','age','black','education','povertyr','bmi','p'))
}
}

