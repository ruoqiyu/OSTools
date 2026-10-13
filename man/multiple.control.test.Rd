\name{multiple.control.test}
\alias{multiple.control.test}
\title{Ordered Comparisons of a Treated Group with Two Control Groups}
\usage{
multiple.control.test(y = NULL, group = NULL, means = NULL,
  se = NULL, n = NULL, treated = NULL, controls = NULL,
  alpha = 0.05, alternative = c("greater", "less"),
  na.rm = FALSE)
}
\arguments{
\item{y, group}{Individual responses and group labels for exactly three independent samples. Alternatively, supply all of \code{means}, \code{se}, and \code{n}. Input ordering and missing-value rules follow \code{contrast.test}.}
\item{means, se, n}{Means, standard errors of the means, and integer sample sizes of at least two. Names align the summary vectors. For unnamed means the labels are \code{group1}, \code{group2}, and \code{group3}.}
\item{treated}{Label identifying the treated group. Required for raw data. For summary input, the default is the first entry of \code{means}. Specifying the label explicitly is recommended.}
\item{controls}{Labels of the two remaining groups, in the order control 1, control 2. The default follows their order in the input.}
\item{alpha}{Level of each reached hypothesis test, strictly between zero and one. The floating-point value \code{1-alpha} must be representable strictly between zero and one. Reject when \code{p.value <= alpha}.}
\item{alternative}{Direction prespecified for the analysis. \code{"greater"} anticipates larger treated outcomes; \code{"less"} reverses all five inequalities.}
\item{na.rm}{Whether to omit raw rows with missing response or group label. Summary inputs must be complete and finite.}
}
\description{Implements the specific five-hypothesis, three-stage testing plan for one treated group and two controls. Uses independent-group Welch contrasts and enforces the stopping rules.}
\details{
Write \eqn{\mu_T}, \eqn{\mu_1}, and \eqn{\mu_2} for the means of the three observed groups. For \code{alternative="greater"}, the nulls are
\deqn{H_1:\mu_T-(\mu_1+\mu_2)/2\leq0,}
\deqn{H_2:\mu_T-\mu_1\leq0,\qquad H_3:\mu_T-\mu_2\leq0,}
\deqn{H_4:\mu_T-2\mu_1+\mu_2\leq0,\qquad H_5:\mu_T+\mu_1-2\mu_2\leq0.}
For \code{"less"}, each null uses the reverse inequality and is tested against a negative contrast. The first contrast always averages the two control means equally, irrespective of their sample sizes; it does not pool control individuals.

Stage 1 tests \eqn{H_1}. Stage 2 is reached only if \eqn{H_1} is rejected; it tests both \eqn{H_2} and \eqn{H_3} and reports either rejection. Stage 3 is reached only if both stage-2 nulls are rejected; it tests \eqn{H_4} and \eqn{H_5}. Rejection of both last nulls supports a treated mean separated from each control mean by more than the controls differ from one another, in the prespecified direction. Later tests are not calculated after a stopping rule applies.

Rosenbaum's ordered procedure uses the logical relationships among these five nulls to control familywise error at \eqn{\alpha} when each component is a valid level-\eqn{\alpha} test. It requires no independence among the component tests. Here the component tests are Welch approximations, so the resulting error control is approximate. This is not a general gatekeeping rule for arbitrary hypotheses or an exact finite-sample guarantee for Welch tests.

The observations in the three samples must be independent. Shared controls, pairing, or clustering require different component variance calculations and are outside this function's scope. The group means are observed comparison-group means; their differences need not be causal effects. Roles, direction, and stopping rules should be specified before inspecting outcomes.
}
\value{A list containing:
\describe{
\item{tests}{Five rows with hypothesis, stage, null inequality, estimate, SE, df, statistic, p.value, log.p.value, tested, rejected, and status. Skipped tests have \code{tested=FALSE}, \code{rejected=NA}, status \code{"skipped"}, and \code{NA} numerical results. Other statuses are \code{"rejected"} and \code{"not rejected"}. P-values are the raw component p-values used in this specific ordered procedure.}
\item{coefficients}{The five-row contrast matrix, with columns named in treated/control order.}
\item{groups, group.summary}{The role labels and input group summaries.}
\item{conclusions}{Logical indicators for the treated comparison with the average, each control, both controls, and the strongest separation conclusion. A skipped comparison supplies no rejection evidence. The sign is given by \code{alternative}.}
\item{stage.reached}{Highest stage actually tested.}
\item{component.tests}{Named \code{htest} results for reached hypotheses only. Their one-sided intervals have confidence level \code{1-alpha}; they are marginal, not simultaneous intervals.}
\item{alpha, alternative, method, fwer.control}{The level, direction, method, and explicit statement that Welch error control is approximate.}
}}
\references{Rosenbaum, P. R. (2008). Testing hypotheses in order. \emph{Biometrika}, 95(1), 248--252. \doi{10.1093/biomet/asm085}.}
\examples{
fit <- multiple.control.test(
  means = c(MW = 14.456, LP = 6.005, SR = 6.400),
  se = c(MW = 0.479, LP = 0.377, SR = 0.143),
  n = c(MW = 67, LP = 57, SR = 61),
  treated = "MW", controls = c("LP", "SR"))
fit$tests[, c("hypothesis", "stage", "p.value", "status")]
fit$conclusions

# A weak first comparison stops the procedure.
multiple.control.test(means = c(T = 0, C1 = 1, C2 = 2),
                      se = c(1, 1, 1), n = c(20, 20, 20))$tests
}
\seealso{\code{\link{contrast.test}}, \code{\link{sens.threesided}}}
