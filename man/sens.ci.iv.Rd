\name{sens.ci.iv}
\alias{sens.ci.iv}
\title{Grid-Based Sensitivity Interval for Paired Instrumental Variables}
\description{Inverts normal approximations to signed-rank sensitivity tests on responses adjusted by a candidate treatment effect.}
\usage{
sens.ci.iv(encouraged.treatment, unencouraged.treatment,
  encouraged.response, unencouraged.response, Gamma=1, alpha=0.05,
  tol=0.01, interval=c(-10,10))
}
\arguments{
\item{encouraged.treatment}{Treatment received by the encouraged member of each pair.}
\item{unencouraged.treatment}{Treatment received by the unencouraged member of each pair.}
\item{encouraged.response}{Responses for the encouraged members.}
\item{unencouraged.response}{Responses for the unencouraged members.}
\item{Gamma}{Sensitivity parameter, at least 1.}
\item{alpha}{One minus the confidence level.}
\item{tol}{Positive grid spacing.}
\item{interval}{Increasing pair of endpoints for the search grid.}
}
\value{A numeric vector of lower and upper endpoints on the specified grid.}
\details{The search is limited to the supplied grid. An endpoint at a search boundary does not establish that the confidence set ends there. With weak instruments the accepted set need not be connected; the returned endpoints summarize its span.}
