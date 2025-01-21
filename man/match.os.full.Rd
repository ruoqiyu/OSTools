\name{match.os.full}
\alias{match.os.full}
\title{Optimal full matching.}
\usage{
match.os.full(z,dist,dat,ncontrol.min=0,ncontrol.max=Inf,exact=NULL,solver='')
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{dist}{ A distance matrix where each treated unit has a row and each control unit has a column.}
\item{dat}{ A data frame with length(z) rows.  If the match is feasible, the matched portion of dat is returned with additional columns that define the match.}
\item{ncontrol.min}{ The minimum ratio of controls to treatments that is to be permitted within a matched set.}
\item{ncontrol.max}{ The maximum ratio of controls to treatments that is to be permitted within a matched set.}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{solver}{ Solver to for the minimum cost network flow problem to find the optimal match. The default solver uses RELAX-IV if rrelaxiv is installed and uses LEMON otherwise.}
}

\value{
A list of results is returned.

\item{data}{ The matched sample, selected rows of dat.}
\item{structure}{ The structure of the match.}
\item{effsize}{ The effective sample size of the match.}
}

\description{
The program finds an optimal full match,  incorporating restrictions on matched sets’ ratios of treatment to control units.
}


\details{
The match minimizes the total distance between treated subjects and their matched controls. 

You MUST install and load the optmatch package to use match.os.full().
}
\examples{
\donttest{
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
head(md)
}
}
\references{
Bertsekas, D. P. and Tseng, P. (1988) The relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.  Fortran and C code: http://www.mit.edu/~dimitrib/home.html.  Available in R via the optmatch package.

Hansen, B. B. (2004). Full matching in an observational study of coaching for the SAT. Journal of the American Statistical Association, 99(467), 609-618.

Hansen, B. B., & Klopfer, S. O. (2006). Optimal full matching and related designs via network flows. Journal of computational and Graphical Statistics, 15(3), 609-627.

}
