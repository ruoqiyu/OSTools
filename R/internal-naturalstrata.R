.natural.require <- function() {
  if (!requireNamespace("natstrat", quietly = TRUE) ||
      utils::packageVersion("natstrat") < "2.0.1")
    stop("Install natstrat (version 2.0.1 or later) to use natural-strata functions.")
}

.natural.backend <- function(expr) {
  # natstrat 2.0.1's caret dummy-variable method looks up contr.ltfr on the
  # search path. Attach its declared dependencies only while making the call.
  before <- search()
  on.exit({
    for (package in setdiff(search(), before))
      detach(package, character.only = TRUE, unload = FALSE)
  }, add = TRUE)
  suppressPackageStartupMessages(base::require("natstrat", character.only = TRUE,
                                               quietly = TRUE))
  force(expr)
}

.natural.inputs <- function(z, stratum, dat, treated) {
  if (!is.data.frame(dat) || !nrow(dat) || !ncol(dat))
    stop("dat must be a nonempty data frame.")
  if (is.null(names(dat)) || anyNA(names(dat)) || any(!nzchar(names(dat))) || anyDuplicated(names(dat)))
    stop("dat must have unique nonempty column names.")
  for (v in list(z = z, stratum = stratum)) {
    if (!is.atomic(v) || !is.null(dim(v)) || length(v) != nrow(dat) || anyNA(v))
      stop("z and stratum must be complete label vectors with one entry per row.")
    if (any(!nzchar(as.character(v)))) stop("Group and stratum labels cannot be empty.")
  }
  z <- if (is.factor(z)) droplevels(z) else factor(z)
  stratum <- if (is.factor(stratum)) droplevels(stratum) else factor(stratum)
  if (nlevels(z) < 2L) stop("At least two groups are required.")
  if (length(treated) != 1L || is.na(treated) || !as.character(treated) %in% levels(z))
    stop("treated must identify one of the group labels.")
  list(z = z, stratum = stratum, treated = as.character(treated),
       counts = table(z, stratum))
}

.natural.flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE.")
}
