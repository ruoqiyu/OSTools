ri.ci <- function(treated.r,control.r,alpha=0.05,alternative='greater',statistic='mean.diff',exact=FALSE,n.mc=1000,range=c(-100,100)) {
  alternative <- match.arg(alternative,c('greater','less','two.sided'))
  statistic <- match.arg(statistic,c('mean.diff','t.test'))
  .ri_validate(treated.r,control.r,statistic,n.mc)
  stopifnot(length(alpha)==1, alpha>0, alpha<1, length(range)==2,
            all(is.finite(range)), range[1]<range[2])
  # Reuse assignments so root finding evaluates a fixed Monte Carlo test.
  assignments <- .ri_assignments(length(treated.r)+length(control.r),length(treated.r),exact,n.mc)
  target <- alpha/if (alternative == 'two.sided') 2 else 1
  testtau <- function(tau, direction) {
    treated <- treated.r-tau
    observed <- .ri_stat(treated,control.r,statistic)
    distribution <- .ri_distribution(c(treated,control.r),assignments,statistic)
    .ri_pvalue(distribution,observed,direction)-target
  }
  lb <- if (alternative == 'less') -Inf else stats::uniroot(testtau,range,direction='greater')$root
  ub <- if (alternative == 'greater') Inf else stats::uniroot(testtau,range,direction='less')$root
  list(lb=lb,ub=ub)
}
