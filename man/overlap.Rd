\name{overlap}
\alias{overlap}
\title{Examine the overlap between the treated and control groups.}
\usage{
overlap(score, z, threshold=0.5)
}
\arguments{
\item{score}{ A vector of scores to examine overlap. It is usually chosen as the logit propensity score or the propensity score.}
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{threshold}{ The multiplier of pooled standard deviation to define overlap.}
}
\value{
\item{no.treated.lack.overlap}{ Number of treated units that lack overlap.}
\item{no.control.lack.overlap}{ Number of control units that lack overlap.}
\item{index.lack.overlap}{ Row indices for the units that lack overlap.}
}
\description{
Examine the overlap between the treated and control groups based on Hansen (2011). Specifically, a treated unit lacks overlap if its logit propensity score is more than threshold (0.5) times a within-group standard deviation above the highest control logit propensity score; a control unit lacks overlap if its logit propensity score is more that threshold (0.5) times a within-group standard deviation below the lowest treated logit propensity score.

}
\examples{
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
X<-cbind(female,age,black,education,povertyr,bmi)
p<-glm(z~female+age+black+education+povertyr+bmi,family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
res=overlap(d$p, d$z)
res$no.treated.lack.overlap
res$no.control.lack.overlap
}
\references{

Hansen, B. B. (2011). Propensity score matching to extract latent experiments from nonexperimental data: A case study. In Looking back: Proceedings of a conference in honor of Paul W. Holland (pp. 149-181). Springer New York.

}
