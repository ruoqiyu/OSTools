match.risk.set <- function(dat,id,time,z,covariates=character(),exact=character(),
                           ncontrol=1,caliper=Inf,eligible=NULL,
                           type=c("initiation","event")) {
  type <- match.arg(type)
  if (!is.data.frame(dat) || !nrow(dat) || anyDuplicated(names(dat)))
    stop("dat must be a nonempty data frame with unique column names.")
  one.name <- function(x) is.character(x) && length(x)==1L && !is.na(x) && x %in% names(dat)
  if (anyDuplicated(c(id,time,z)) || !all(vapply(list(id,time,z),one.name,logical(1)))) stop("id, time, and z must name columns in dat.")
  for (v in list(covariates,exact)) if (!is.character(v) || anyNA(v) || anyDuplicated(v) || any(!v %in% names(dat)))
    stop("covariates and exact must be vectors of unique column names.")
  if (any(c(id,time,z) %in% c(covariates,exact))) stop("Design columns cannot also be covariates or exact variables.")
  reserved <- c("mset",".risk.time",".risk.row")
  if (any(reserved %in% names(dat))) stop("dat contains reserved matching-output columns.")
  ids <- as.character(dat[[id]]);times <- dat[[time]];zz <- dat[[z]]
  if (anyNA(ids) || any(!nzchar(ids))) stop("id must contain nonempty, nonmissing identifiers.")
  if (!(is.numeric(times) || inherits(times,"Date")) || any(!is.finite(as.numeric(times))))
    stop("time must be finite numeric decision times or Date values.")
  if (!(is.numeric(zz) || is.logical(zz)) || anyNA(zz) || any(!zz %in% c(0,1)))
    stop("z must identify the event at each decision with zero or one.")
  if (anyDuplicated(dat[c(id,time)])) stop("There must be at most one row per person and decision time.")
  if (length(ncontrol)!=1L || !is.numeric(ncontrol) || !is.finite(ncontrol) || ncontrol<1 || ncontrol!=floor(ncontrol))
    stop("ncontrol must be a positive integer.")
  if (length(caliper)!=1L || !is.numeric(caliper) || is.na(caliper) || caliper<0)
    stop("caliper must be nonnegative (Inf is allowed).")
  ok <- rep(TRUE,nrow(dat))
  if (!is.null(eligible)) {
    if (!one.name(eligible)) stop("eligible must be NULL or a column name.")
    e <- dat[[eligible]]
    if (!(is.logical(e) || is.numeric(e)) || anyNA(e) || any(!e %in% c(0,1)))
      stop("The eligibility column must contain nonmissing logical or 0-1 values.")
    ok <- as.logical(e)
  }
  if (anyNA(dat[c(covariates,exact)])) stop("Matching covariates and exact variables must not be missing.")
  for (nm in covariates) if (is.numeric(dat[[nm]]) && any(!is.finite(dat[[nm]])))
    stop("Numeric matching covariates must be finite.")
  if (type=="initiation" && any(table(ids[zz==1])>1L))
    stop("In initiation mode each person can initiate only once; z must mark initiation, not ongoing treatment.")
  used <- character();past <- character();pieces <- list();diagnostics <- list()
  reason <- rep(NA_character_,nrow(dat));reason[zz==1 & !ok] <- "ineligible"
  mset <- 0L
  decision.times <- sort(unique(times))
  for (kk in seq_along(decision.times)) {
    tt <- decision.times[kk]
    at <- which(times==tt)
    pool <- at[ok[at] & !ids[at] %in% used & (type=="event" | !ids[at] %in% past)]
    # Canonical ordering makes input-row permutations immaterial, including ties.
    pool <- pool[order(ids[pool],method="radix")]
    ti <- pool[zz[pool]==1];co <- pool[zz[pool]==0]
    reason[at[zz[at]==1 & ok[at] & ids[at] %in% used]] <- "already matched"
    if (length(ti)) reason[ti] <- "insufficient eligible controls within caliper"
    done <- 0L
    if (length(ti) && length(co)) {
      # Estimate distance only from the current eligible risk set. Neither
      # future rows nor unused outcome columns enter this calculation.
      xx <- .os_covariates(dat[pool,covariates,drop=FALSE])
      if (ncol(xx)) {
        xx <- apply(xx,2,rank);xx <- matrix(xx,nrow=length(pool))
        cv <- stats::cov(xx)
        rat <- sqrt(stats::var(seq_along(pool))/diag(cv))
        cv <- diag(rat,nrow=length(rat)) %*% cv %*% diag(rat,nrow=length(rat))
        inv <- MASS::ginv(cv)
        d <- t(vapply(match(ti,pool),function(k)
          sqrt(pmax(0,stats::mahalanobis(xx[match(co,pool),,drop=FALSE],xx[k,],inv,inverted=TRUE))),numeric(length(co))))
      } else d <- matrix(0,length(ti),length(co))
      dim(d) <- c(length(ti),length(co))
      for (nm in exact) d[outer(as.character(dat[[nm]][ti]),as.character(dat[[nm]][co]),"!=")] <- Inf
      d[d>caliper] <- Inf
      reason[ti[rowSums(is.finite(d))>=ncontrol]] <- "not selected in maximum-cardinality match"
      # Solve each exact stratum independently. Distances retain the covariance
      # estimated from the whole current eligible pool, as in earlier releases.
      cell <- if (length(exact)) as.integer(interaction(dat[pool,exact,drop=FALSE],
        drop=TRUE,lex.order=TRUE)) else rep.int(1L,length(pool))
      assignments <- list()
      for (cc in sort(unique(cell))) {
        it <- which(cell[match(ti,pool)]==cc)
        ic <- which(cell[match(co,pool)]==cc)
        if (!length(it) || length(ic)<ncontrol) next
        a <- .os_risk_assignment(d[it,ic,drop=FALSE],ncontrol)
        for (j in which(lengths(a)>0L)) assignments[[as.character(it[j])]] <- ic[a[[j]]]
      }
      for (j in sort(as.integer(names(assignments)))) {
        choice <- assignments[[as.character(j)]]
        selected <- c(ti[j],co[choice]);mset <- mset+1L;done <- done+1L
        part <- dat[selected,,drop=FALSE];part$mset <- mset
        part$.risk.time <- times[selected];part$.risk.row <- selected
        pieces[[length(pieces)+1L]] <- part
        used <- c(used,ids[selected]);reason[ti[j]] <- NA_character_
      }
    }
    diagnostics[[length(diagnostics)+1L]] <- data.frame(time=tt,eligible=length(pool),
      treated=length(ti),controls=length(co),sets=done,unmatched=length(ti)-done)
    past <- union(past,ids[at[zz[at]==1]])
  }
  md <- if (length(pieces)) do.call(rbind,pieces) else {
    a <- dat[FALSE,,drop=FALSE];a$mset <- integer();a$.risk.time <- times[FALSE];a$.risk.row <- integer();a
  }
  rownames(md) <- NULL
  selected <- seq_len(nrow(dat)) %in% md$.risk.row
  u <- which(zz==1 & !selected)
  unmatched <- data.frame(row=u,id=ids[u],time=times[u],reason=reason[u],stringsAsFactors=FALSE)
  ans <- list(data=md,z=zz[md$.risk.row],selected=selected,
      unmatched.treated=unmatched,diagnostics=do.call(rbind,diagnostics),
      ncontrol=ncontrol,type=type,algorithm="Chronological optimal rank-Mahalanobis matching without person reuse",
      call=match.call())
  class(ans) <- "os_risk_set_match";ans
}

# Binary assignment with all-or-none treated groups. First maximize the number
# of complete 1:k sets, then minimize total edge distance at that cardinality.
# A single penalized objective or cloned treated slots need not preserve this
# lexicographic priority and can produce incomplete groups under hard calipers.
.os_risk_assignment <- function(d,k) {
  nt <- nrow(d);nc <- ncol(d)
  ans <- rep(list(integer()),nt)
  if (!nt || nc<k) return(ans)
  edge <- which(is.finite(d) & rowSums(is.finite(d))>=k,arr.ind=TRUE)
  ne <- nrow(edge)
  if (!ne) return(ans)
  # x_tc selects an allowed edge; b_t selects a complete treated group.
  # sum_c x_tc = k*b_t; sum_t x_tc <= 1 for every control.
  mat <- slam::simple_triplet_matrix(
    i=c(edge[,1],nt+edge[,2],seq_len(nt)),
    j=c(seq_len(ne),seq_len(ne),ne+seq_len(nt)),
    v=c(rep.int(1,2L*ne),rep.int(-k,nt)),
    nrow=nt+nc,ncol=ne+nt)
  dirs <- c(rep.int("==",nt),rep.int("<=",nc))
  rhs <- c(rep.int(0,nt),rep.int(1,nc))
  solve <- function(obj,mat,dirs,rhs,max=FALSE) {
    fit <- Rglpk::Rglpk_solve_LP(obj,mat,dirs,rhs,
      types=rep.int("B",ne+nt),max=max,control=list(presolve=TRUE))
    if (fit$status!=0L || any(!is.finite(fit$solution)))
      stop("Optimal risk-set assignment failed; no greedy fallback was used.")
    fit$solution
  }
  first <- solve(c(rep.int(0,ne),rep.int(1,nt)),mat,dirs,rhs,max=TRUE)
  count <- sum(round(first[ne+seq_len(nt)]))
  if (!count) return(ans)
  mat <- slam::simple_triplet_matrix(
    i=c(mat$i,rep.int(nt+nc+1L,nt)),j=c(mat$j,ne+seq_len(nt)),
    v=c(mat$v,rep.int(1,nt)),nrow=nt+nc+1L,ncol=ne+nt)
  best <- solve(c(d[edge],rep.int(0,nt)),mat,c(dirs,"=="),c(rhs,count))
  chosen <- edge[best[seq_len(ne)]>.5,,drop=FALSE]
  sizes <- tabulate(chosen[,1],nbins=nt)
  if (any(sizes!=0 & sizes!=k) || sum(sizes==k)!=count ||
      anyDuplicated(chosen[,2]) || any(abs(best-round(best))>1e-6))
    stop("Optimal risk-set assignment returned an invalid solution.")
  for (j in which(sizes>0L)) ans[[j]] <- sort(chosen[chosen[,1]==j,2])
  ans
}

print.os_risk_set_match <- function(x,...) {
  cat(x$algorithm,"\n")
  cat("Mode:",x$type,"; matched sets:",length(unique(x$data$mset)),
      "; matched people:",nrow(x$data),"; unmatched treated rows:",nrow(x$unmatched.treated),"\n")
  invisible(x)
}
