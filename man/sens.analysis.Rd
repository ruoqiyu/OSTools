\name{sens.analysis}
\alias{sens.analysis}
\title{
Sensitivity analysis in observational studies.
}
\usage{
sens.analysis(y=NULL,D=NULL,Tobs=NULL,Gamma=1,tau=0,alternative='two-sided',
method='m',inner=0,trim=3,lambda=1/2,weight.par=c(1,1,1))
}
\arguments{
\item{y}{
An n by J matrix where (i) the rows are n matched sets, (ii) the first column is the treated response in a set, columns 2 to J contain the responses of controls in the same matched set. Every set must have J-1 controls, and NAs are not allowed in y.
}
\item{D}{
An input argument required when method='mcnemar': the number of discordant pairs.
}
\item{Tobs}{
Observed value of test statistic for the McNemar test: number of discordant pairs in which treated unit has a 1.
}
\item{Gamma}{
The sensitivity parameter: Gamma=1 for a randomization test, Gamma>1 for sensitivity bounds. This parameter Gamma is denoted by the upper case Greek letter Gamma in the literature, e.g., Rosenbaum (2007, 2014).
}
\item{tau}{
If tau=0, test the null hypothesis of no treatment effect.  If tau is not 0, test the null hypothesis that the treatment effect is an additive shift of tau against the alternative that the effect is larger than tau in the direction specified by alternative.
}
\item{alternative}{
If alternative='greater', the one-sided alternative hypothesis is that treatment increases the level of response by at least tau. If alternative='less', the one-sided alternative that the treatment decreases the level of response by at least tau. If alternative='two-sided', the alternative hypothesis is two-sided.
}
\item{method}{
If method='m', using weighted Huber-Maritz M-statistics determined by the parameters, namely inner, trim, lambda, weight.par. For example:

(i) For Huber's unweighted m-statistic, set inner=0, trim=3, lambda=1/2, weight.par=c(1,1,1). The unweighted method is often a good choice in small samples when there are a few pairs or sets (say 20 sets) or when the number of controls in each matched set is 6 or more.

(ii) For a sturdy and all-purpose weighed M-statistic, set inner=0, trim=3, lambda=1/2, weight.par=c(20,12,20). This weighted method often works better than the unweighted version with 2-4 controls per matched set.  

(iii) For matched pairs, set inner = 1/2, trim =2, lambda = 1/2, weight.par=(1,1,1). This method performs well for matched pairs, as seen in the evaluations in Rosenbaum (2013).

(iv) For short tailed distributions, such as the Normal distribution, set inner=0, trim=3, lambda=1/2, weight.par=c(20,16,20). This method often performs better for short-tailed Normal errors than the all-purpose weights.

(v) For long tailed distributions, such as the t-distribution with 4 degrees of freedom, set inner=0, trim=3, lambda=1/2, weight.par=c(20,12,19). These weights redescend. This method often performs better for long-tailed Normal errors than the all-purpose weights.

(vi) The Quade ranks sets using ordinary ranks (1, 2, ..., n) applied to ranges of M-scores within sets, in parallel with Quade (1979) and Tardiff (1987).  It sets inner=0, trim=3, lambda=1/2, weight.par=c(2,2,2)

(vii) The permutational t-test is unweighted and permutes the observations themselves without ranking or scoring. It sets inner=0, trim=Inf, lambda=1/2, weight.par=c(1,1,1).   This method is the permutation test that uses the treated-minus-control difference in means as a test statistic.

if method='mcnemar', using the McNemar test for paired binary outcomes.
}
\item{inner}{
Inner trimming to increase design sensitivity. Use of inner<0 or inner>trim will generate an error.
}
\item{trim}{
Outer trimming for resistance to outliers. If trim = Inf, there is no trimming.
}
\item{lambda}{
Observations are scaled by the lambda quantile of the absolute pair differences, defaulting to the median of all paired absolute differences; see Rosenbaum (2007, section 4.2) for a precise definition in the case of multiple controls. 
}
\item{weight.par}{
A vector of three parameters (m,m1,m2) that determine the weights, as in expression (5) in Rosenbaum (2014) or expression (8) in Rosenbaum (2011) where they are called m (for m), underline(m) for m1, and overline(m) for m2.
}

}
\value{
\item{pval}{Approximate upper bound on the P-value. }
\item{deviate}{Deviate that is compared to the upper tail of the standard Normal distribution to obtain the P-value.}
\item{statistic}{Value of the test statistic.}
\item{expectation}{Maximum null expectation of the test statistic for the given value of Gamma.}
\item{variance}{Among null distributions that yield the maximum expectation, variance is the maximum possible variance for the given value of Gamma.}
}

\description{
Computes the large sample approximation to the upper bound on the P-value testing the null hypothesis of no treatment effect versus the alternative hypothesis in a matched observational study with a fixed number of controls matched to each treated subject. If method='m', the function uses an unweighted or a weighted version of Huber-Maritz M-statistics as test statistics, with weights defined by weights.par.
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
o<-match.os(d$z, dist, d, fine=d$education, ncontrol=1)
md<-o$data
sens.analysis(cbind(md$homocysteine[md$z==1],md$homocysteine[md$z==0]),Gamma=1)
sens.analysis(cbind(md$homocysteine[md$z==1],md$homocysteine[md$z==0]),Gamma=1.5)
}
}

\references{

Huber, P. (1981) Robust Statistics.  New York: Wiley, 1981.

Maritz, J. S. (1979) Exact robust confidence intervals for location.  Biometrika 1979, 66, 163-166. 

Gastwirth, J. L., Krieger, A. M., and Rosenbaum, P. R. (2000) Asymptotic separability in sensitivity analysis.  Journal of the Royal Statistical Society B 2000, 62, 545-556.

Quade, D. (1979). Using weighted rankings in the analysis of complete blocks with additive block effects. Journal of the American Statistical Association, 74(367), 680-683.

Rosenbaum, P. R. (2007) Sensitivity analysis for m-estimates, tests and confidence intervals in matched observational studies.  Biometrics, 2007, 63, 456-464. 

Rosenbaum, P. R. (2013) Impact of multiple matched controls on design sensitivity in observational studies.  Biometrics, 2013, 69, 118-127.

Rosenbaum, P. R. (2014) Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls.  Journal of the American Statistical Association, 109(507), 1145-1158

Rosenbaum, P. R. (2015). Two R packages for sensitivity analysis in observational studies. Observational Studies, 1(2), 1-17.  The Observational Studies journal is available free on-line.

Rosenbaum, P. R. (2018). Sensitivity analysis for stratified comparisons in an observational study of the effect of smoking on homocysteine
levels. The Annals of Applied Statistics, 12(4), 2312-2334.

Tardif, S. (1987). Efficiency and optimality results for tests based on weighted rankings. Journal of the American Statistical Association, 82(398), 637-644.

}




