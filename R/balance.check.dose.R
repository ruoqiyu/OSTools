balance.check.dose <- function(X, dose, mset) {
  # A single, deterministic split rule is used before and after matching.
  split_dose <- function(rows) {
    rows <- rows[order(dose[rows],rows)]
    cut <- floor(length(rows)/2L)
    list(lower=rows[seq_len(cut)],higher=rows[seq.int(cut+1L,length(rows))])
  }
  group_means <- function(rows) {
    groups <- split_dose(rows)
    cbind(higher=colMeans(X[groups$higher,,drop=FALSE]),
          lower=colMeans(X[groups$lower,,drop=FALSE]))
  }
  standardize <- function(difference, scale) {
    result <- difference/scale
    result[scale==0 & difference==0] <- 0
    result
  }

  if (is.numeric(X) && is.null(dim(X))) X <- matrix(X,ncol=1L)
  if (!(is.matrix(X) || is.data.frame(X)))
    stop("X must be a numeric matrix, data frame, or vector.",call.=FALSE)
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X)<4L || ncol(X)<1L || any(!is.finite(X)))
    stop("X must contain at least four rows and one column of finite numeric covariates.",call.=FALSE)
  n <- nrow(X)
  if (!is.numeric(dose) || !is.null(dim(dose)) || length(dose)!=n || any(!is.finite(dose)))
    stop("dose must be a finite numeric vector with one value per row of X.",call.=FALSE)
  if ((!is.atomic(mset) && !is.factor(mset)) || !is.null(dim(mset)) || length(mset)!=n)
    stop("mset must be a vector with one set label per row of X.",call.=FALSE)
  keep <- !is.na(mset)
  if (!any(keep)) stop("mset contains no matched subjects.",call.=FALSE)
  sets <- split(which(keep),as.character(mset[keep]),drop=TRUE)
  if (any(lengths(sets)<2L)) stop("Each matched set must contain at least two subjects.",call.=FALSE)
  if (any(vapply(sets,function(rows) length(unique(dose[rows]))==1L,logical(1))))
    warning("Some matched sets have no dose variation; tied doses are split by original row order.",call.=FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("X",seq_len(ncol(X)))

  before.groups <- split_dose(seq_len(n))
  before <- group_means(seq_len(n))
  scale <- sqrt((apply(X[before.groups$higher,,drop=FALSE],2L,stats::var)+
                 apply(X[before.groups$lower,,drop=FALSE],2L,stats::var))/2)
  # Equal set weights follow the course helper's averaging convention.
  after <- Reduce(`+`,lapply(sets,group_means))/length(sets)
  result <- cbind(before[,"higher"],before[,"lower"],
                  after[,"higher"],after[,"lower"],
                  standardize(before[,"higher"]-before[,"lower"],scale),
                  standardize(after[,"higher"]-after[,"lower"],scale))
  colnames(result) <- c("Before Higher Mean","Before Lower Mean",
                        "After Higher Mean","After Lower Mean","Before SMD","After SMD")
  rownames(result) <- colnames(X)
  result
}
