.ri_stat <- function(treated, control, statistic) {
  difference <- mean(treated)-mean(control)
  if (statistic == 'mean.diff') return(difference)
  m <- length(treated)
  n <- length(control)
  variance <- ((m-1)*stats::var(treated)+(n-1)*stats::var(control))/(m+n-2)
  se <- sqrt(variance*(1/m+1/n))
  if (se == 0 && difference == 0) return(0)
  difference/se
}

.ri_validate <- function(treated, control, statistic, n.mc) {
  stopifnot(is.numeric(treated), is.numeric(control), length(treated)>0,
            length(control)>0, all(is.finite(treated)), all(is.finite(control)))
  if (statistic == 't.test') stopifnot(length(treated)>1, length(control)>1)
  stopifnot(length(n.mc)==1, is.finite(n.mc), n.mc>=1, n.mc==floor(n.mc))
}

.ri_assignments <- function(n, m, exact, n.mc) {
  if (exact) return(utils::combn(n,m))
  matrix(replicate(n.mc,sample.int(n,m)),nrow=m)
}

.ri_distribution <- function(r, assignments, statistic) {
  vapply(seq_len(ncol(assignments)),function(i) {
    treated <- assignments[,i]
    .ri_stat(r[treated],r[-treated],statistic)
  },numeric(1))
}

.ri_pvalue <- function(distribution, observed, alternative) {
  # A union of both tails avoids counting zero statistics twice.
  switch(alternative,
         greater=mean(distribution>=observed),
         less=mean(distribution<=observed),
         two.sided=mean(abs(distribution)>=abs(observed)))
}

ri.test <- function(treated.r,control.r,alternative='greater',statistic='mean.diff',exact=FALSE,n.mc=1000) {
  alternative <- match.arg(alternative,c('greater','less','two.sided'))
  statistic <- match.arg(statistic,c('mean.diff','t.test'))
  .ri_validate(treated.r,control.r,statistic,n.mc)
  r <- c(treated.r,control.r)
  observed <- .ri_stat(treated.r,control.r,statistic)
  assignments <- .ri_assignments(length(r),length(treated.r),exact,n.mc)
  distribution <- .ri_distribution(r,assignments,statistic)
  pval <- .ri_pvalue(distribution,observed,alternative)
  if (exact) return(list(pval=pval,teststat=distribution,obs.test.stat=observed))
  error <- 1.96*sqrt(pval*(1-pval)/n.mc)
  list(pval=pval,lowerci=max(0,pval-error),upperci=min(1,pval+error),
       mc.test.stat=distribution,obs.test.stat=observed)
}
