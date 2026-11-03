\name{partial.conjunction}
\alias{partial.conjunction}
\alias{print.os_partial_conjunction}
\title{Partial-Conjunction Tests and Prespecified Individual Follow-Up}
\usage{
partial.conjunction(p, method = c("fisher", "truncated"),
  trunc = 0.2, alpha = 0.05, r.followup = NULL, log.p = FALSE)
\method{print}{os_partial_conjunction}(x, digits = 4, ...)
}
\arguments{
\item{p}{Numeric vector of component p-values, or their natural logarithms when \code{log.p = TRUE}. The input rules are those of \code{combine.pvalues}. Unique names label the individual hypotheses.}
\item{method}{Fisher or truncated-product combination.}
\item{trunc}{Prespecified cutoff in \eqn{(0,1]}. Fisher's method uses 1.}
\item{alpha}{Test level strictly between zero and one.}
\item{r.followup}{Optional integer from 1 to \code{length(p)}, chosen before examining these p-values, specifying the gate for individual follow-up. \code{NULL} reports count claims only. Do not choose this target using the returned count or the most favorable p-value.}
\item{log.p}{Whether the inputs are log probabilities.}
\item{x}{A returned \code{os_partial_conjunction} object.}
\item{digits}{Number of digits used when printing results.}
\item{...}{Additional arguments passed to \code{print.data.frame}.}
}
\description{Tests hypotheses about how many component nulls are false, with a fixed-sequence adjustment for count claims and optional gated individual tests at a prespecified target.}
\details{
For \eqn{n} component hypotheses, the partial-conjunction null for \eqn{r} states that at most \eqn{r-1} component nulls are false. Order the p-values from smallest to largest and apply \code{combine.pvalues} to the \eqn{n-r+1} largest, for each \eqn{r=1,\ldots,n}. The returned raw p-values need not be increasing in \eqn{r}.

To report a lower bound on the number of false nulls, test \eqn{r=1,2,\ldots,n} and stop at the first nonrejection. The adjusted p-value at \eqn{r} is the cumulative maximum of the raw values through \eqn{r}. Thus the count is not the largest isolated \eqn{r} whose raw p-value happens to be small. All comparisons are computed using log probabilities; \code{p.adjusted} here denotes this fixed-sequence adjustment, not Holm adjustment.

If \code{r.followup} is supplied and is less than \eqn{n}, first test that prespecified partial-conjunction null using its raw p-value at level \code{alpha}. Only if this gate rejects are individual p-values compared with \eqn{\alpha/(n-r)}. When \eqn{r=n}, rejection of the gate supports all component nulls being false and no division or further individual threshold is needed. The gate is distinct from the optional sequence of count claims. A failed gate yields no individual rejections. A count of \eqn{r} effects need not identify \eqn{r} specific effects.

The component null p-values must be independent, or satisfy the required evidence-factor joint stochastic bounds under every relevant intersection null. This is stronger than marginal validity and is not checked by the function. The combination and truncation cutoff are fixed in advance. The follow-up target must also be prespecified; setting it to the estimated count does not justify the reported familywise error guarantee. See Karmakar and Small (2020), Proposition 5, for the gated rule.

The count-only fixed sequence and the prespecified gate with its individual follow-up are two separate level-\code{alpha} reporting plans. The latter controls errors in its one prespecified count claim and the individual tests; it does not automatically control the union with every count claim from the separate fixed sequence. Choose the reporting plan in advance, or provide additional multiplicity control if both families of claims are to be reported together. Both are returned to make the distinction inspectable.

The function analyzes one vector at a time. For a sensitivity grid, call it separately on each row of a matrix of valid component bounds. These are pointwise claims at the supplied bias bounds. Simultaneous inference over arbitrary parameter searches is not provided; structured-range results require additional assumptions and a suitable testing plan.
}
\value{An object of class \code{os_partial_conjunction}. Its \code{tests} data frame contains \code{r}, raw \code{p.value}, fixed-sequence \code{p.adjusted}, sequence \code{reject}, \code{log.p.value}, and \code{log.p.adjusted}. \code{lower.bound} is the number supported by the ordered sequence. Input probabilities, their logs, the sorted hypothesis labels, \code{alpha}, \code{combination}, \code{trunc}, and \code{assumptions} are retained. \code{followup} is \code{NULL} unless requested; otherwise it contains the prespecified \code{r}, raw gate p-value and log p-value, \code{gate.reject}, \code{threshold} and \code{log.threshold} (both \code{NA} for \eqn{r=n}), and an input-order \code{tests} data frame of hypothesis labels, ordinary/log p-values and individual \code{reject} decisions.}
\references{
Benjamini, Y. and Heller, R. (2008). Screening for partial conjunction hypotheses. \emph{Biometrics}, 64, 1215--1222.

Wang, J. and Owen, A. B. (2019). Admissibility in partial conjunction testing. \emph{Journal of the American Statistical Association}, 114, 158--168. \doi{10.1080/01621459.2017.1385465}.

Karmakar, B. and Small, D. S. (2020). Assessment of the extent of corroboration of an elaborate theory of a causal hypothesis using partial conjunctions of evidence factors. \emph{The Annals of Statistics}, 48, 3283--3311. \doi{10.1214/19-AOS1929}.
}
\examples{
# Chapter 13: five lead-study comparisons.
p <- c(treated = 0.00007, exposure = 0.004, dose = 0.10,
       hygiene = 0.009, hygiene2 = 0.42)
ans <- partial.conjunction(p, r.followup = 3)
ans$tests
ans$lower.bound
ans$followup$tests

# A count can be supported without identifying that many individual effects.
p2 <- c(0.004007, 0.025994, 0.252909, 0.031257, 0.611224)
ans2 <- partial.conjunction(p2, method = "truncated", r.followup = 2)
ans2$lower.bound
ans2$followup$tests
}
\seealso{\code{\link{combine.pvalues}}}
