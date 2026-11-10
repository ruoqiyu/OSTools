.tobit.inputs <- function(y,z,mset,Gamma,method,inner,trim,lambda) {
  scalar <- function(v) is.numeric(v) && length(v)==1L && is.null(dim(v)) && is.finite(v)
  if (!is.numeric(y) || !is.null(dim(y)) || !length(y) ||
      any(!is.finite(y)) || any(y<0 | y>1)) stop("y must be a finite numeric vector in [0,1].")
  if (!(is.numeric(z) || is.logical(z)) || !is.null(dim(z)) ||
      length(z)!=length(y) || anyNA(z) || any(!z %in% c(0,1)))
    stop("z must be a 0-1 vector with the same length as y.")
  if (!is.atomic(mset) || !is.null(dim(mset)) || length(mset)!=length(y) || anyNA(mset))
    stop("mset must contain one nonmissing label per response.")
  if (!scalar(Gamma) || Gamma<1) stop("Gamma must be a finite scalar at least one.")
  if (!scalar(inner) || inner<0 || !is.numeric(trim) || length(trim)!=1L ||
      !is.null(dim(trim)) || is.na(trim) || trim<=0 || trim<inner ||
      !scalar(lambda) || lambda<=0 || lambda>=1)
    stop("Require 0 <= inner <= trim, trim > 0, and 0 < lambda < 1.")
  groups <- split(seq_along(y),factor(mset,levels=unique(mset)))
  if (any(vapply(groups,function(i) length(i)<2L || sum(z[i])!=1,logical(1))))
    stop("Each matched set must contain exactly one treated individual and at least one control.")
  groups <- lapply(groups,function(i) c(i[z[i]==1],i[z[i]==0]))
  list(y=y,z=as.integer(z),mset=mset,groups=groups,Gamma=Gamma,
       method=method,inner=inner,trim=trim,lambda=lambda)
}

.tobit.scores <- function(input,tau) {
  a <- pmax(0,input$y-(1-input$z)*tau)
  ys <- lapply(input$groups,function(i) a[i])
  diffs <- lapply(ys,function(v) outer(v,v,"-"))
  scale <- 1; fallback <- FALSE
  robust <- input$method=="m" && (is.finite(input$trim) || input$inner>0)
  if (robust) {
    ad <- unlist(lapply(diffs,function(d) abs(d[row(d)!=col(d)])),use.names=FALSE)
    scale <- as.numeric(stats::quantile(ad,input$lambda,names=FALSE))
    if (scale==0) {
      pos <- ad[ad>0]
      if (length(pos)) {scale <- min(pos);fallback <- TRUE} else scale <- 1
    }
  }
  q <- lapply(diffs,function(d) {
    if (robust) {
      w <- abs(d)/scale
      magnitude <- if (!is.finite(input$trim)) pmax(0,w-input$inner) else
        if (input$trim==input$inner) as.numeric(w>input$inner) else
          pmin(1,pmax(0,w-input$inner)/(input$trim-input$inner))
      dim(magnitude) <- dim(d)
      d <- sign(d)*magnitude
    }
    rowMeans(d)
  })
  qm <- matrix(NA_real_,length(q),max(lengths(q)))
  for (i in seq_along(q)) qm[i,seq_along(q[[i]])] <- q[[i]]
  list(scores=qm,adjusted=a,scale=scale,scale.fallback=fallback)
}

.tobit.tail <- function(q,Gamma) {
  # The backend supplies separable worst-case moments. Sum the Normal tail
  # directly, rather than subtracting its CDF from one.
  fit <- sensitivityfull::separable1f(q,gamma=Gamma)
  v <- fit$variance
  if (v<0 && abs(v)<1e-12*max(1,sum(q^2,na.rm=TRUE))) v <- 0
  if (!is.finite(v) || v<0) stop("Could not compute nonnegative sensitivity variance.")
  gap <- fit$statistic-fit$expectation
  if (v==0) {
    dev <- if (abs(gap)<=1e-12*max(1,abs(fit$statistic),abs(fit$expectation))) 0 else sign(gap)*Inf
    lp <- if (dev>0) -Inf else 0
  } else {
    dev <- gap/sqrt(v)
    lp <- stats::pnorm(dev,lower.tail=FALSE,log.p=TRUE)
  }
  list(p.value=exp(lp),log.p.value=lp,deviate=dev,statistic=fit$statistic,
       expectation=fit$expectation,variance=v,score=gap)
}

.tobit.evaluate <- function(input,tau,alternative) {
  sc <- .tobit.scores(input,tau)
  # A greater reduction means a smaller treated outcome.
  greater <- .tobit.tail(-sc$scores,input$Gamma)
  less <- .tobit.tail(sc$scores,input$Gamma)
  ans <- if (alternative=="greater") greater else if (alternative=="less") less else {
    a <- if (greater$log.p.value<=less$log.p.value) greater else less
    a$log.p.value <- min(0,log(2)+min(greater$log.p.value,less$log.p.value))
    a$p.value <- exp(a$log.p.value);a
  }
  ans$pval <- ans$p.value
  ans$tau <- tau;ans$Gamma <- input$Gamma;ans$alternative <- alternative
  ans$method <- input$method;ans$adjusted <- sc$adjusted
  ans$scale <- sc$scale;ans$scale.fallback <- sc$scale.fallback
  ans$scores <- c(greater=greater$score,less=less$score)
  ans$approximation <- "Separable Normal approximation"
  ans
}

.tobit.regions <- function(evaluate,grid,columns,tol,epsilon) {
  vals <- vapply(grid,evaluate,numeric(3))
  roots <- numeric()
  for (j in columns) for (i in seq_len(length(grid)-1L)) {
    a <- vals[j,i];b <- vals[j,i+1L]
    if (is.finite(a) && is.finite(b) && sign(a)*sign(b)<0) {
      rr <- stats::uniroot(function(t) evaluate(t)[j],c(grid[i],grid[i+1L]),
                           tol=min(tol/10,1e-10))$root
      roots <- c(roots,rr)
    }
  }
  nodes <- sort(unique(c(grid,roots)))
  # Treat explicitly scanned intervals independently; do not presume a single
  # monotone crossing. This remains a numerical scan, not an exhaustive proof.
  keep <- function(t) all(evaluate(t)[columns]>= -epsilon)
  pieces <- list()
  for (i in seq_len(length(nodes)-1L)) {
    if (keep((nodes[i]+nodes[i+1L])/2)) pieces[[length(pieces)+1L]] <- nodes[c(i,i+1L)]
  }
  for (t in nodes) if (keep(t)) pieces[[length(pieces)+1L]] <- c(t,t)
  if (!length(pieces)) return(matrix(numeric(),0,2,dimnames=list(NULL,c("lower","upper"))))
  pieces <- do.call(rbind,pieces);pieces <- pieces[order(pieces[,1],pieces[,2]),,drop=FALSE]
  merged <- list(pieces[1,])
  if (nrow(pieces)>1L) for (i in 2:nrow(pieces)) {
    k <- length(merged)
    if (pieces[i,1]<=merged[[k]][2]+tol) merged[[k]][2] <- max(merged[[k]][2],pieces[i,2])
    else merged[[k+1L]] <- pieces[i,]
  }
  result <- do.call(rbind,merged);colnames(result) <- c("lower","upper");result
}
