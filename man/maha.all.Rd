\name{maha.all}
\alias{maha.all}
\title{Creates a robust Mahalanobis distance among all units.}
\usage{
maha.all(X)
}
\arguments{
\item{X}{ A matrix with length(z) rows giving the covariates.  X should be of full column rank.}
}

\value{
A distance matrix for each pair of units.
}
\description{
Computes a robust Mahalanobis distance matrix for use in non-bipartite matching and crossmatch test for covariate balance evaluation. In this case, we compute the distance for all possible pairs of two units (treated and treated, treated and control, control and control).

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
dist<-maha.all(X)
dist[1:4,1:6]
detach(nh0506Homocysteine)

}
\references{

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

}
