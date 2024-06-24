\name{ri.ci}
\alias{ri.ci}
\title{Confidence interval for randomization inference}
\usage{
ri.ci(treated.r,control.r,alpha=0.05,alternative='greater',
statistic='mean.diff',exact=FALSE,n.mc=1000,range=c(-100,100))
}
\arguments{
\item{treated.r}{ A vector of responses in the treated group.}
\item{control.r}{ A vector of responses in the control group.}
\item{alpha}{ Coverage probability of confidence interval is 1-alpha.}
\item{alternative}{ Direction of the alternative hypothesis. It can be one of 'greater', 'less', and 'two.sided'.}
\item{statistic}{ Choice of test statistic. It can be one of 'mean.diff' and 't.test'.}
\item{exact}{ If exact=TRUE, calculate the exact p-value. If exact=FALSE, calculate the p-value based on monte-carlo simulations}
\item{n.mc}{ Number of monte-carlo simulations to approximate the true randomization distribution of the test statistic.}
\item{range}{ Range to search for the root in the uniroot function, which is used to determine the two end points of the confidence interval.}
}
\value{
\item{lb}{Lower bound of the (1-alpha) confidence interval.}
\item{ub}{Upper bound of the (1-alpha) confidence interval.}
}
\description{
Confidence interval for the treatment effect, which is constructed by inverting  the randomization tests.
}
\examples{
treated.r.fishoil=c(8,12,10,14,2,0,0);
control.r.fishoil=c(-6,0,1,2,-3,-4,2);

ri.ci(treated.r=treated.r.fishoil,control.r=control.r.fishoil,alpha=0.05,alternative='two.sided',
      statistic='mean.diff',exact=T)
}

