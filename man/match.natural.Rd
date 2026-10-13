\name{match.natural}
\alias{match.natural}
\title{Construct a Natural-Strata Comparison Sample}
\usage{
match.natural(z, stratum, dat, balance.formulas, ratio = 1,
  treated = 1, denom.variance = c("pooled", "treated"),
  weight.by.size = 0, seed = NULL, runs = 1,
  solver = c("Rglpk", "gurobi"), correct.sizes = TRUE, ...)
}
\arguments{
\item{z}{Complete group labels, one per row of \code{dat}. Factors retain their observed level order; other labels are sorted.}
\item{stratum}{Complete exact-stratum labels, one per row. Every stratum must contain the treated group and enough observations of each control group.}
\item{dat}{Data frame in the original observation order. A working copy gains the factor column \code{stratum}; an existing column of that name must agree with the supplied labels. Do not include outcomes in balance formulas.}
\item{balance.formulas}{A formula or list of formulas in the syntax of \code{natstrat::generate_constraints}. Covariates occur on the left; \code{1} requests overall balance and \code{stratum} requests within-stratum balance on the right. Categorical covariates are expanded into indicators.}
\item{ratio}{Positive controls per treated observation. A scalar applies to every control group. For different ratios, supply a uniquely named vector for all controls or all groups; the treated group's ratio is one. Names align the ratios with group labels.}
\item{treated}{The treated group label, default \code{1}.}
\item{denom.variance}{Standardization using the mean of original group variances or the original treated variance. The backend falls back to pooled variance when treated variance is zero.}
\item{weight.by.size}{Between zero and one. Zero uses the backend's default count scaling; one gives greater relative importance to smaller strata.}
\item{seed}{Optional nonnegative integer seed for randomized rounding.}
\item{runs}{Positive integer number of backend rounding runs.}
\item{solver}{Backend solver. \code{"Rglpk"} is the open-source default; \code{"gurobi"} requires a separately installed licensed solver.}
\item{correct.sizes}{Whether to request the backend's correction of realized counts. Realized counts and ratios must still be inspected.}
\item{...}{Named additional arguments to \code{natstrat::optimize_controls}, such as \code{time_limit} or supplementary \code{treated_star} and \code{q_star_s}. Primary design arguments and \code{q_s} cannot be replaced here. Supplementary results remain in \code{backend}.}
}
\description{An OSTools interface to the natural-strata optimization and randomized rounding in the optional natstrat package (version 2.0.1 or later). Retains the treated group and selects controls to balance pretreatment covariates.}
\details{
Target counts are \code{round(ratio * treated count)} in each stratum, using R's rounding convention. Infeasible targets and zero target counts raise errors rather than silently discarding treated observations. Integer rounding may prevent exact proportionality. \code{proportional} checks whether realized group ratios are common across strata; \code{target.ratio.achieved} separately checks whether all requested ratios were attained. Only exact common realized ratios justify replacing treated-count-weighted stratum differences by an unweighted difference of group means.

The backend handles missing covariates as documented by \code{generate_constraints}; decisions about missing data should precede design construction. Selection does not establish ignorability. Inspect overall and within-stratum balance, and use appropriate stratified inference. A fixed seed improves reproducibility within the same backend/solver environment; equivalent optimization solutions can differ across environments. Temporary backend dependency attachments are removed on return, including errors.
}
\value{A list containing \code{data} (selected rows in original order), \code{selected} (logical vector in original row order), selected \code{z} and \code{stratum}, realized \code{counts} and \code{ratios}, named \code{target.ratios}, rounded \code{target.counts}, logical \code{proportional} and \code{target.ratio.achieved}, the \code{treated} label, generated \code{constraints}, the original optimizer result in \code{backend}, \code{seed}, and \code{solver}. Count and ratio tables have group rows and stratum columns.}
\examples{
if (requireNamespace("natstrat", quietly = TRUE)) {
  z <- rep(c(1, 0), c(4, 8))
  st <- rep(c("A", "B"), c(2, 2))[c(1:4, 1:4, 1:4)]
  dat <- data.frame(x = c(1, 3, 2, 4, 0, 2, 1, 3, 4, 6, 5, 7))
  fit <- match.natural(z, st, dat, x ~ 1 + stratum,
                       ratio = 1, seed = 7)
  fit$counts
  fit$proportional
  balance.check.strat(z, dat, st, selected = fit$selected)
}
}
\seealso{\code{\link{balance.check.strat}}, \code{\link[natstrat]{generate_constraints}}, \code{\link[natstrat]{optimize_controls}}}
