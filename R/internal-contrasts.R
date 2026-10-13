.contrast.labels <- function(labels, what) {
  if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))
    stop(what, " must have nonmissing, nonempty, unique group names.")
  labels
}

.contrast.align <- function(x, labels, what) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) != length(labels) ||
      any(!is.finite(x)))
    stop(what, " must be a finite numeric vector with one entry per group.")
  if (!is.null(names(x))) {
    .contrast.labels(names(x), what)
    if (!setequal(names(x), labels))
      stop(what, " names must match the group names exactly.")
    x <- x[labels]
  }
  stats::setNames(as.numeric(x), labels)
}

.contrast.groups <- function(y, group, means, se, n, na.rm) {
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm))
    stop("na.rm must be TRUE or FALSE.")
  raw <- !is.null(y) || !is.null(group)
  if (raw) {
    if (!is.null(means) || !is.null(se) || !is.null(n))
      stop("Supply raw y and group, or means, se and n; do not mix them.")
    if (is.null(y) || is.null(group)) stop("Supply both y and group.")
    if (!is.numeric(y) || !is.null(dim(y)) || length(y) == 0L)
      stop("y must be a numeric response vector.")
    if (!is.atomic(group) || !is.null(dim(group)) || length(group) != length(y))
      stop("group must be a vector of labels with the same length as y.")
    missing <- is.na(y) | is.na(group)
    if (any(missing) && !na.rm)
      stop("Missing responses or group labels require na.rm=TRUE.")
    y <- y[!missing]
    group <- group[!missing]
    if (!length(y) || any(!is.finite(y)))
      stop("Retained responses must be finite and nonempty.")
    labels <- if (is.factor(group)) levels(droplevels(group)) else unique(as.character(group))
    .contrast.labels(labels, "group")
    group <- as.character(group)
    n <- vapply(labels, function(g) sum(group == g), numeric(1))
    if (any(n < 2)) stop("Each group must contain at least two observations.")
    means <- vapply(labels, function(g) mean(y[group == g]), numeric(1))
    se <- vapply(labels, function(g) {
      x <- y[group == g]
      centered <- x - means[g]
      if (all(is.finite(centered))) {
        scale <- max(abs(centered))
        if (scale == 0) return(0)
        scale * (stats::sd(centered / scale) / sqrt(n[g]))
      } else {
        # Subtraction itself may overflow for opposite, very large responses.
        scale <- max(abs(x))
        scale * (stats::sd(x / scale) / sqrt(n[g]))
      }
    }, numeric(1))
  } else {
    if (is.null(means) || is.null(se) || is.null(n))
      stop("Supply y and group, or all of means, se and n.")
    if (!is.numeric(means) || !is.null(dim(means)) || !length(means) ||
        any(!is.finite(means)))
      stop("means must be a nonempty finite numeric vector.")
    labels <- names(means)
    if (is.null(labels)) {
      if (!is.null(names(se)) || !is.null(names(n)))
        stop("Name means to align named se and n vectors.")
      labels <- paste0("group", seq_along(means))
    }
    .contrast.labels(labels, "means")
  }
  means <- .contrast.align(means, labels, "means")
  se <- .contrast.align(se, labels, "se")
  n <- .contrast.align(n, labels, "n")
  if (any(se < 0)) stop("se contains standard errors of group means and must be nonnegative.")
  if (any(n < 2 | n != floor(n))) stop("n must contain integer sample sizes of at least two.")
  list(means = means, se = se, n = n, labels = labels, raw = raw,
       summary = data.frame(group = labels, mean = unname(means),
                            se = unname(se), n = unname(n), row.names = NULL))
}
