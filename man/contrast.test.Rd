\name{contrast.test}
\alias{contrast.test}
\title{Welch Test and Interval for an Independent-Group Contrast}
\usage{
contrast.test(y = NULL, group = NULL, contrast, means = NULL,
  se = NULL, n = NULL, null = 0,
  alternative = c("two-sided", "greater", "less"),
  conf.level = 0.95, na.rm = FALSE)
}
\arguments{
\item{y}{Numeric vector of individual responses, for independent samples. Supply either \code{y} and \code{group}, or all of \code{means}, \code{se}, and \code{n}.}
\item{group}{Group labels of the same length as \code{y}. Factor levels determine the group order; otherwise the order of first appearance is used. Unused factor levels are dropped.}
\item{contrast}{Finite numeric coefficients, one per group, with at least one nonzero coefficient. Named coefficients are matched to group labels. Unnamed coefficients follow the order of \code{means}, or the raw-data group order. Naming coefficients is recommended. A general linear combination is allowed; coefficients need not sum to zero.}
\item{means}{Finite numeric group means. Unique names identify groups. Unnamed means receive labels \code{group1}, \code{group2}, and so on.}
\item{se}{Nonnegative standard errors of the group means, not standard deviations of individual responses. Named vectors are aligned to \code{means}; unnamed vectors follow its order. Named \code{se} or \code{n} requires named \code{means}.}
\item{n}{Integer group sample sizes of at least two. Each raw-data group must also have at least two retained observations.}
\item{null}{Finite null value for the linear combination of group means.}
\item{alternative}{Two-sided or directional alternative. The R spelling \code{"two.sided"} is also accepted.}
\item{conf.level}{Confidence level strictly between zero and one. One-sided alternatives return the corresponding one-sided interval.}
\item{na.rm}{Whether to omit raw-data rows with missing response or group label. Infinite responses are not permitted. Summary inputs must always be complete and finite.}
}
\description{Tests a linear combination of means from independent groups using a Welch--Satterthwaite approximation. Accepts individual responses or group summaries.}
\details{
Let \eqn{a_g} be the coefficients and \eqn{e_g} the standard errors of group means. The estimated contrast is \eqn{\sum_g a_g\bar Y_g}, its variance estimate is \eqn{V=\sum_g a_g^2 e_g^2}, and the degrees of freedom are
\deqn{\nu=V^2/\sum_g\{(a_g^2 e_g^2)^2/(n_g-1)\}.}
Scaling the variance contributions before squaring avoids unnecessary numerical overflow and underflow. Constant groups are permitted when the overall contrast standard error is positive. A contrast with zero standard error cannot be tested and raises an error.

P-values use the appropriate Student t tail directly on the log scale. \code{log.p.value} remains informative when a finite tail probability is too small to represent as an ordinary double. Two-sided intervals use two-sided quantiles; directional intervals have one infinite endpoint. Finite interval endpoints are checked, and unrepresentable intervals raise an error suggesting rescaling.

The independence assumption excludes paired observations, overlapping control samples, clustered samples, and correlated group summaries unless the required covariance has been handled outside this function. These intervals are approximate and marginal; they are not simultaneous intervals for multiple contrasts. Differences between observational group means need additional assumptions to be interpreted causally.
}
\value{An object of class \code{htest}, containing \code{statistic}, \code{parameter} (degrees of freedom), \code{p.value}, \code{log.p.value}, \code{conf.int}, \code{estimate}, \code{null.value}, \code{stderr}, \code{alternative}, \code{method}, and \code{data.name}. It also includes named \code{coefficients} and a \code{group.summary} data frame. The returned \code{alternative} uses standard \code{htest} spelling \code{"two.sided"}.}
\examples{
# Lead-zinc miners minus the equal-weight average of the controls.
means <- c(MW = 14.456, LP = 6.005, SR = 6.400)
se <- c(MW = 0.479, LP = 0.377, SR = 0.143)
n <- c(MW = 67, LP = 57, SR = 61)
contrast.test(means = means, se = se, n = n,
              contrast = c(MW = 1, LP = -0.5, SR = -0.5),
              alternative = "greater")

# Raw independent samples; named coefficients avoid ordering ambiguity.
y <- c(1, 2, 4, 5, 2, 3, 6, 8, 9)
group <- c(rep("treated", 4), rep("control", 5))
contrast.test(y, group, contrast = c(treated = 1, control = -1))
}
\seealso{\code{\link{multiple.control.test}}, \code{\link[stats]{t.test}}}
