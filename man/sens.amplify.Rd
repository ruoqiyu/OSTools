\name{sens.amplify}
\alias{sens.amplify}
\title{Interpret a Sensitivity Parameter by Amplification}
\usage{sens.amplify(Gamma, Lambda = NULL)}
\arguments{
\item{Gamma}{A nonempty numeric vector of finite sensitivity parameters at least 1.}
\item{Lambda}{Optional finite treatment-association odds ratios. For Gamma greater than 1, each Lambda must exceed its corresponding Gamma. For Gamma equal to 1, Lambda must be at least 1. Supply equal-length vectors or one scalar. If omitted, use the symmetric point Lambda = Delta.}
}
\description{Returns amplification pairs and the within-pair treatment probability bounds used to interpret a sensitivity analysis.}
\details{For Gamma greater than 1, the curve satisfies Gamma = (Lambda * Delta + 1)/(Lambda + Delta). The function computes Delta = (Gamma * Lambda - 1)/(Lambda - Gamma). The symmetric point is Gamma + sqrt(Gamma squared - 1). For Gamma = 1, a supplied Lambda gives Delta = 1; the default representative is (1,1). Other points on the degenerate Gamma = 1 curve have Lambda = 1 and Delta greater than 1. Finite points are returned; infinite asymptotes are not included. Amplification originates in matched-pair sensitivity analysis, and does not generally represent an arbitrary confounding model for every design.}
\value{A data frame with columns \code{Gamma}, \code{Lambda}, \code{Delta}, \code{probability.lower}, \code{probability.upper}, and \code{symmetric}. The probability limits are 1/(1+Gamma) and Gamma/(1+Gamma). The symmetric column gives the common association at Lambda = Delta.}
\examples{
sens.amplify(c(1, 1.5, 2, 5))
sens.amplify(5, Lambda = c(6, 8, 10, 20))
sens.amplify(c(1, 2, 5), Lambda = 20)
}
\references{Rosenbaum, P. R. and Silber, J. H. (2009). Amplification of sensitivity analysis in observational studies. Journal of the American Statistical Association, 104, 1398--1405. \doi{10.1198/jasa.2009.tm08470}.}
\seealso{\code{\link{sens.value}}, \code{\link{sens.analysis}}}
