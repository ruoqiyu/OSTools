\name{balance.check.full}
\alias{balance.check.full}
\title{Check standardized mean differences (SMDs) of full matching.}
\usage{
balance.check.full(fdata,mdata,fz,mz,mset,fmissing=0,mmissing=0,parameter="ate")
}
\arguments{
\item{fdata}{ A full data frame with length(fz) rows and columns being variables that need to check SMDs. fdata and mdata must have the same variables with same column names in the same order.}
\item{mdata}{ A matched data frame with length(mz) rows and columns being variables that need to check SMDs. fdata and mdata must have the same variables with same column names in the same order.}
\item{fz}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control for subjects in the full data set.}
\item{mz}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control for subjects in the matched data set.}
\item{mset}{ A vector of matched set information.}
\item{fmissing}{ A matrix of missingness indicators for the full data.}
\item{mmissing}{ A matrix of missingness indicators for the matched data.}
\item{parameter}{ It can be "ate" or "att", which corresponds to different weights when evaluating the balance after matching.}
}
\value{ The output is a list of three parts.
\item{var}{ Names of covariates for which the balance is evaluated.}
\item{smd.before}{ The standardized mean differences before matching.}
\item{smd.after}{ The standardized mean differences after matching.}
}
\description{
The function is used to create a table of mean and SMDs to check the balance before and after full matching.
}
\examples{
# To run this example, you must load the optmatch package.
# The optmatch is available on CRAN and Github.
\donttest{
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
X<-cbind(female,age,black,education,povertyr,bmi)
p<-glm(z~female+age+black+education+povertyr+bmi,family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
dist<-maha.dense(d$z,X,matrix=TRUE)
o<-match.os.full(d$z, dist, d)
matcheddata<-o$data
Xm<-subset(matcheddata, select=c('female','age','black','education','povertyr','bmi','p'))
balance.check.full(cbind(X,p),Xm,d$z,matcheddata$z,matcheddata$mset)
}
}
\references{

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

}
