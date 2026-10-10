match.dose <- function(X, dose, method=c("full", "pair"),
                       distance=c("rank", "mahalanobis", "euclidean"),
                       lambda=0, near.exact=NULL, eta=1, mindist=0.01,
                       penalty=100000, ndiscard=0) {
  # All matching-specific helpers are private to this function.
  check_number <- function(value, name, allow.inf=FALSE) {
    if (!is.numeric(value) || length(value)!=1L || is.na(value) ||
        value<0 || (!allow.inf && !is.finite(value)))
      stop(name, " must be one nonnegative number",
           if (!allow.inf) " (finite)" else "", ".", call.=FALSE)
  }
  construct_dist_mat <- function(x, metric) {
    if (metric == "rank") return(maha.all(x))
    if (metric == "euclidean") return(as.matrix(stats::dist(x)))
    varying <- apply(x, 2L, function(v) length(unique(v))>1L)
    x <- x[, varying, drop=FALSE]
    if (!ncol(x)) return(matrix(0, nrow(x), nrow(x)))
    inverse <- MASS::ginv(stats::cov(x))
    result <- vapply(seq_len(nrow(x)), function(i)
      stats::mahalanobis(x, x[i, ], inverse, inverted=TRUE), numeric(nrow(x)))
    pmax((result+t(result))/2, 0)
  }
  solve_pairs <- function(cost) {
    ans <- nbpMatching::nonbimatch(nbpMatching::distancematrix(cost), precision=8)$halves
    edges <- cbind(as.integer(ans$Group1.Row), as.integer(ans$Group2.Row))
    if (!identical(sort(as.integer(edges)), seq_len(nrow(cost))) ||
        any(!is.finite(cost[edges])))
      stop("No valid match was returned. Relax the matching restrictions.", call.=FALSE)
    edges
  }
  # Recover an edge cover from either layer of the doubled graph. Removing
  # redundant edges leaves disjoint stars, including pairs, even with zero ties.
  extract_match_set <- function(pairs, cost, nearest, offset) {
    within <- pairs[,1L]>offset & pairs[,1L]<=offset+n &
              pairs[,2L]>offset & pairs[,2L]<=offset+n
    edges <- pairs[within, , drop=FALSE]-offset
    unmatched <- pairs[pairs[,1L]<=n & pairs[,2L]>n, 1L]
    if (length(unmatched)) edges <- rbind(edges, cbind(unmatched, nearest[unmatched]))
    edges <- unique(cbind(pmin(edges[,1L],edges[,2L]),pmax(edges[,1L],edges[,2L])))
    degree <- tabulate(as.integer(edges), nbins=n)
    for (k in order(cost[edges], decreasing=TRUE)) {
      if (all(degree[edges[k,]]>1L)) {
        degree[edges[k,]] <- degree[edges[k,]]-1L
        edges[k,] <- 0L
      }
    }
    edges <- edges[edges[,1L]>0L, , drop=FALSE]
    if (any(tabulate(as.integer(edges), nbins=n)==0L))
      stop("Internal error: the edge cover omitted a subject.", call.=FALSE)
    edges
  }
  label_sets <- function(edges) {
    # Union-find makes labels independent of the solver's pair order.
    parent <- seq_len(n)
    root <- function(i) {
      while (parent[i]!=i) i <- parent[i]
      i
    }
    for (k in seq_len(nrow(edges))) {
      a <- root(edges[k,1L]); b <- root(edges[k,2L])
      parent[max(a,b)] <- min(a,b)
    }
    used <- sort(unique(as.integer(edges)))
    result <- rep(NA_integer_,n)
    roots <- vapply(used,root,integer(1))
    result[used] <- match(roots,sort(unique(roots)))
    result
  }
  summarize_sets <- function(edges, mset, raw.distance, costs) {
    ids <- sort(unique(mset[!is.na(mset)]))
    do.call(rbind,lapply(ids,function(id) {
      subjects <- which(mset==id)
      e <- edges[mset[edges[,1L]]==id,,drop=FALSE]
      degree <- tabulate(as.integer(e),nbins=n)
      center <- subjects[which.max(degree[subjects])]
      if (max(degree[subjects])!=length(subjects)-1L)
        stop("Internal error: a recovered matched set is not a star.", call.=FALSE)
      subdist <- raw.distance[subjects,subjects,drop=FALSE]
      gaps <- abs(dose[e[,1L]]-dose[e[,2L]])
      data.frame(mset=id, size=length(subjects), center=center,
                 mean.dist=mean(subdist[upper.tri(subdist)]),
                 star.dist=mean(raw.distance[e]), dose.diff=mean(gaps),
                 min.dose.diff=min(gaps), n.close=sum(gaps<mindist),
                 cost=sum(costs[e]))
    }))
  }

  method <- match.arg(method)
  distance <- match.arg(distance)
  if (is.numeric(X) && is.null(dim(X))) X <- matrix(X,ncol=1L)
  if (!(is.matrix(X) || is.data.frame(X)))
    stop("X must be a numeric matrix, data frame, or vector.",call.=FALSE)
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X)<2L || ncol(X)<1L || any(!is.finite(X)))
    stop("X must contain at least two rows and one column of finite numeric covariates.",call.=FALSE)
  n <- nrow(X)
  if (!is.numeric(dose) || !is.null(dim(dose)) || length(dose)!=n || any(!is.finite(dose)))
    stop("dose must be a finite numeric vector with one value per row of X.",call.=FALSE)
  check_number(lambda,"lambda"); check_number(eta,"eta")
  check_number(mindist,"mindist"); check_number(penalty,"penalty",TRUE)
  check_number(ndiscard,"ndiscard")
  if (ndiscard!=floor(ndiscard) || ndiscard>n-2L)
    stop("ndiscard must be an integer between zero and nrow(X)-2.",call.=FALSE)
  if (method=="full" && ndiscard!=0)
    stop("ndiscard is only available for pair matching.",call.=FALSE)
  if (method=="pair" && lambda!=0)
    stop("lambda is only available for full matching.",call.=FALSE)
  if (method=="pair" && (n-ndiscard)%%2L!=0L)
    stop("Pair matching requires an even number of retained subjects; adjust ndiscard.",call.=FALSE)
  if (!is.null(near.exact)) {
    if (is.character(near.exact)) {
      if (is.null(colnames(X)) || anyNA(match(near.exact,colnames(X))))
        stop("near.exact must name columns of X.",call.=FALSE)
      near.exact <- match(near.exact,colnames(X))
    }
    if (!is.numeric(near.exact) || anyNA(near.exact) ||
        any(near.exact!=floor(near.exact) | near.exact<1 | near.exact>ncol(X)))
      stop("near.exact must contain column names or valid column numbers.",call.=FALSE)
    near.exact <- unique(near.exact)
  }

  raw.distance <- construct_dist_mat(X,distance)
  diag(raw.distance) <- 0
  cost <- raw.distance
  near.penalty <- eta*max(raw.distance)
  for (j in near.exact) {
    mismatch <- outer(X[,j],X[,j],"!=")
    cost[mismatch] <- cost[mismatch]+near.penalty
  }
  too.close <- abs(outer(dose,dose,"-"))<mindist
  cost[too.close] <- cost[too.close]+penalty
  diag(cost) <- 0
  if (any(is.nan(cost)) || any(cost<0))
    stop("The distance calculation overflowed; rescale the covariates or penalties.",call.=FALSE)

  if (method=="pair") {
    augmented <- matrix(0,n+ndiscard,n+ndiscard)
    augmented[seq_len(n),seq_len(n)] <- cost
    if (ndiscard>0L) {
      phantom <- n+seq_len(ndiscard)
      augmented[phantom,phantom] <- Inf
      diag(augmented) <- 0
    }
    pairs <- solve_pairs(augmented)
    edges <- pairs[pairs[,1L]<=n & pairs[,2L]<=n,,drop=FALSE]
  } else {
    available <- cost
    diag(available) <- Inf
    nearest <- apply(available,1L,which.min)
    minimum <- available[cbind(seq_len(n),nearest)]
    if (any(!is.finite(minimum)))
      stop("At least one subject has no admissible partner; relax mindist or penalty.",call.=FALSE)
    bridge <- matrix(Inf,n,n)
    diag(bridge) <- 2*minimum+2*lambda
    doubled <- rbind(cbind(cost,bridge),cbind(bridge,cost))
    pairs <- solve_pairs(doubled)
    if (any(pairs[pairs[,1L]<=n & pairs[,2L]>n,2L]-n !=
            pairs[pairs[,1L]<=n & pairs[,2L]>n,1L]))
      stop("The solver selected a forbidden cross-layer connection.",call.=FALSE)
    first <- extract_match_set(pairs,cost,nearest,0L)
    second <- extract_match_set(pairs,cost,nearest,n)
    objective <- function(e) sum(cost[e])+lambda*(2*nrow(e)-n)
    edges <- if (objective(first)<=objective(second)) first else second
  }
  edges <- cbind(pmin(edges[,1L],edges[,2L]),pmax(edges[,1L],edges[,2L]))
  edges <- edges[order(edges[,1L],edges[,2L]),,drop=FALSE]
  mset <- label_sets(edges)
  diagnostics <- summarize_sets(edges,mset,raw.distance,cost)
  selected <- data.frame(i=edges[,1L],j=edges[,2L],mset=mset[edges[,1L]],
                         distance=raw.distance[edges],cost=cost[edges],
                         dose.diff=abs(dose[edges[,1L]]-dose[edges[,2L]]))
  if (any(diagnostics$n.close>0L))
    warning("Some selected connections have dose differences below mindist; the finite penalty is a soft restriction.",call.=FALSE)
  list(mset=mset, dose=dose, edges=selected, diagnostics=diagnostics,
       mean.distance.before=mean(raw.distance[upper.tri(raw.distance)]),
       discarded=which(is.na(mset)), method=method, distance=distance,
       objective=sum(selected$cost)+if (method=="full") lambda*sum(diagnostics$size-2L) else 0,
       call=match.call())
}
