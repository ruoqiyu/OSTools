\name{ri.test}
\alias{ri.test}
\title{Randomization inference for testing no treatment effect}
\usage{
ri.test(treated.r,control.r,alternative='greater',statistic='mean.diff',
exact=FALSE,n.mc=1000)
}
\arguments{
\item{treated.r}{ A vector of responses in the treated group.}
\item{control.r}{ A vector of responses in the control group.}
\item{alternative}{ Direction of the alternative hypothesis. It can be one of 'greater', 'less', and 'two.sided'.}
\item{statistic}{ Choice of test statistic. It can be one of 'mean.diff' and 't.test'.}
\item{exact}{ If exact=TRUE, calculate the exact p-value. If exact=FALSE, calculate the p-value based on monte-carlo simulations}
\item{n.mc}{ Number of monte-carlo simulations to approximate the true randomization distribution of the test statistic.}
}
\value{
If exact=TRUE, the output is a list of p-value, all possible test statistic values under the randomization distribution, and the observed test statistic value. If exact=FALSE, the output is a list of p-value, its confidence interval, all test statistic values under the monte-carlo simulations, and the observed test statistic value. 
}
\description{
Randomization test for no treatment effect using the difference in sample means between the treated and control subjects or the two-sample t test as the test statistic.
}
\examples{
treated.r.fishoil=c(8,12,10,14,2,0,0);
control.r.fishoil=c(-6,0,1,2,-3,-4,2);

# Test no treatment effect for fish oil data
res=ri.test(treated.r.fishoil,control.r.fishoil,
            alternative='greater',statistic='mean.diff',exact=T)
res$pval

set.seed(10)
res=ri.test(treated.r.fishoil,control.r.fishoil,alternative='greater',
            statistic='mean.diff',exact=F,n.mc=1000)
res$pval
res$lowerci
res$upperci
}

