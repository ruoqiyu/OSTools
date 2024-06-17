\name{maha.sparse}
\alias{maha.sparse}
\title{Creates a robust Mahalanobis distance for matching based on a sparse network.}
\usage{
maha.sparse(z,X,p=rep(1,length(z)),caliper=1,stdev=FALSE,constant=NULL,ncontrol=1,
exact=NULL,nearexact=NULL,penalty=100,subX=NULL,ties.all=TRUE)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. }
\item{X}{ A matrix with length(z) rows giving the covariates.  X should be of full column rank.}
\item{p}{ A vector of length(z)=length(p) giving the variable used to define the caliper.  Typically, p is the propensity score or its rank.}
\item{caliper}{ If two individuals differ on p by more than caliper, we will not calculate the distance for this pair. If caliper is a positive number, then a symmetric caliper is applied. If caliper is a vector of a negative number and a positive number, then an asymmetric caliper is applied.}
\item{stdev}{ If stdev = TRUE, caliper is interpreted in units of an equally weighted pooled
standard deviation; that is,caliper is replaced by caliper*sp where sp is
sqrt((var(dx[z==1])+var(dx[z==0]))/2).}
\item{constant}{ If the number of pairs within a caliper is greater than constant, we will select the constant closest ones.}
\item{ncontrol}{ A positive integer giving the number of controls to be matched to each treated subject.  If ncontrol is too large, the match will be infeasible.}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{nearexact}{ If not NULL, then a vector of length length(z) giving variable that need to be exactly matched. }
\item{penalty}{ The penalty for a mismatch on nearexact.}
\item{subX}{ If a subset matching is required, the variable that the subset matching is based on. That is, for each level of subX, extra treated will be discarded in order to have the number of matched treated subjects being the minimum size of treated and control groups. If exact matching on a variable x is desired and discarding extra treated is fine if there are more treated than controls for a certain level k, set exact=x, subX=x.}
\item{ties.all}{ If ties.all is True, include all ties while choosing nearest neighbors. In this case, some treated may have more than constant controls. Otherwise, randomly select one or several controls to make sure there are not more than constant controls for each treated.}
}

\value{
\item{d}{ A distance list for each pair within the caliper distance and constant constraint.}
\item{start}{ The treated subject for each distance.}
\item{end}{ The control subject for each distance.}
}
\description{
Computes a robust Mahalanobis distance list for use in sparse matching. In this case, we will only calculate the distance for pairs within the caliper on p. If the caliper is too small, the matching may be infeasible. For the smallest caliper that keeps feasiblibity, refer to optcal() in package 'bigmatch'.

This function and its use are discussed in Rosenbaum (2010). It is preferred when the dataset is large. The robust
Mahalanobis distance in described in Chapter 8 of Rosenbaum (2010).

}
\details{
The usual Mahalanobis distance works well for multivariate Normal covariates, but can exhibit odd behavior with typical covariates. Long tails or an outlier in a covariate can yield a large estimated
variance, so the usual Mahalanobis distance pays little attention to large differences in this covariate. Rare binary covariates have a small variance, so a mismatch on a rare binary covariate is
viewed by the usual Mahalanobis distance as extremely important. If you were matching for binary covariates indicating US state of residence, the usual Mahalanobis distance would regard a
mismatch for Wyoming as much worse than a mismatch for California.

The robust Mahalanobis distance uses ranks of covariates rather than the covariates themselves, but the variances of the ranks are not adjusted for ties, so ties do not make a variable more important. Binary covariates are, of course, heavily tied.
}
\examples{
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
X<-cbind(female,age,black,education,povertyr,bmi)
p<-glm(z~female+age+black+education+povertyr+bmi,family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
\donttest{
#apply symmetric caliper 0.15 on propensity score
dist1<-maha.sparse(d$z,X,p,0.15)
length(dist1$d)
}
#apply asymmetric caliper c(-0.2,0.1) on propensity score
dist2<-maha.sparse(d$z,X,p,c(-0.2,0.1))
length(dist2$d)
}
\references{

Yu, R., Silber, J. H., & Rosenbaum, P. R. (2020).  Matching methods for observational studiesderived from large administrative databases (with Discussion).Statistical Science, 35(3), 338-355.

}
