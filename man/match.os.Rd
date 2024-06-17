\name{match.os}
\alias{match.os}
\title{Minimum-distance near-fine matching.}
\usage{
match.os(z,dist,dat,p=rep(1,length(z)),exact=NULL,fine=rep(1,length(z)),ncontrol=1,
penalty=ifelse(is.matrix(dist),round(max(dist)*1000),round(max(dist$d)*1000)),
s.cost=100,subX=NULL,solver='rlemon')
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{dist}{ A distance object with three components: d, start, end, typically created by maha_dense or maha_sparse. d[i] gives the distance between the (start[i])th treated and the (end[i]-sum(z))th control.}
\item{dat}{ A data frame with length(z) rows.  If the match is feasible, the matched portion of dat is returned with additional columns that define the match.}
\item{p}{ A vector of length(z)=length(p) giving the variable used to define the caliper.  Typically, p is the propensity score or its rank. If the dense match is performed, use the default p=rep(1,length(z)).}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{fine}{ A vector of with length(z)=length(fine) giving the nominal levels that are to be nearly-finely balanced.}
\item{ncontrol}{ A positive integer giving the number of controls to be matched to each treated subject.  If ncontrol is too large, the match will be infeasible.}
\item{penalty}{ A numeric penalty imposed for each violation of fine balance.}
\item{s.cost}{ The scaling factor for cost of the each pair of treated and control while rounding the cost.}
\item{subX}{ If a subset matching is required, the variable that the subset matching is based on. That is, for each level of subX, extra treated will be discarded in order to have the number of matched treated subjects being the minimum size of treated and control groups. If exact matching on a variable x is desired and discarding extra treated is fine if there are more treated than controls for a certain level k, set exact=x, subX=x.}
\item{solver}{ Solver to for the minimum cost network flow problem to find the optimal match. The default solver is rlemon. An alternative solver is rrelaxiv.}
}

\value{
If the match is infeasible, a warning is issued.  Otherwise, a list of results is returned.

A match may be infeasible if the caliper or constant on p is too small, or ncontrol is too large, or if exact matching for exact is impossible.

\item{feasible}{ Indicator of whether matching is feasible or not.}
\item{data}{ The matched sample, selected rows of dat.}
\item{x}{ A vector of indicators of whether each treated-control pair is included in the matched sample.}
}

\description{
The program finds an optimal near-fine match with a given caliper on p, plus directional penalties on dx to offset the distribution imbalances. That is, it finds a near-fine match that minimizes the penalized Mahalanobis distance. In order to avoid the distortion of the original distribution by large penalties, it has the option of apply asymmetric calipers on those covariates.
}


\details{
The match minimizes the total distance between treated subjects and their matched controls subject to a near-fine balance constraint imposed as a penalty on imbalances. Another set of directional penalties on dx can be imposed in order to offset the distribution imbalances. In order to avoid the case of matching far pairs to get close means, the user can the option of apply asymmetric calipers on covariates in dx. We add a larger penalty for pairs outside the asymmetric caliper to avoid infeasibility issue. But a match may be infeasible if the caliper on p is too small. In this case, increase the caliper, or find the smallest caliper that gives a feasible matching by using optcal() in package 'bigmatch'.

For discussion of networks for fine-balance, see Rosenbaum (1989, Section 3) and Rosenbaum (2010).
For near-fine balannce balance, see Yang et al. (2012).

You MUST install and load the optmatch package to use match().
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
dist<-add.caliper(dist, d$z, d$p, c(-.5,.15), stdev=TRUE, penalty=1000)
o<-match.os(d$z, dist, d, fine=d$education, ncontrol=2)
md<-o$data
head(md)
}
}
\references{
Bertsekas, D. P. and Tseng, P. (1988) The relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.  Fortran and C code: http://www.mit.edu/~dimitrib/home.html.  Available in R via the optmatch package.

Rosenbaum, P.R. (1989) Optimal matching in observational studies. Journal of the American Statistical Association, 84, 1024-1032.

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

Yang, D., Small, D. S., Silber, J. H., and Rosenbaum, P. R. (2012) Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes.  Biometrics, 68, 628-636.

Yu, R., & Rosenbaum, P. R. (2019). Directional penalties for optimal matching in observationalstudies.Biometrics, 75(4), 1380-1390.
}
