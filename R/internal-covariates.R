# Internal: numeric design matrix, with mean imputation and missingness indicators.
.os_covariates <- function(x) {
  if (is.null(dim(x))) x <- data.frame(x=x)
  x <- as.data.frame(x)
  n <- nrow(x)
  columns <- list()
  for (i in seq_along(x)) {
    v <- x[[i]]
    if (is.character(v) || is.factor(v)) {
      v <- factor(v)
      missing <- is.na(v)
      for (level in levels(v)) {
        columns[[length(columns)+1L]] <- as.numeric(!missing & v == level)
      }
      if (any(missing)) columns[[length(columns)+1L]] <- as.numeric(missing)
    } else {
      stopifnot(is.numeric(v) || is.logical(v))
      missing <- is.na(v)
      if (any(missing)) {
        columns[[length(columns)+1L]] <- as.numeric(missing)
        v[missing] <- if (all(missing)) 0 else mean(v[!missing])
      }
      stopifnot(all(is.finite(v)))
      columns[[length(columns)+1L]] <- as.numeric(v)
    }
  }
  if (!length(columns)) return(matrix(numeric(0), n, 0))
  result <- do.call(cbind, columns)
  varying <- apply(result,2,function(v) length(unique(v)) > 1L)
  result[,varying,drop=FALSE]
}
