match.natural <- function(z, stratum, dat, balance.formulas, ratio = 1,
                          treated = 1, denom.variance = c("pooled", "treated"),
                          weight.by.size = 0, seed = NULL, runs = 1,
                          solver = c("Rglpk", "gurobi"), correct.sizes = TRUE, ...) {
  .natural.require()
  input <- .natural.inputs(z, stratum, dat, treated)
  z <- input$z; stratum <- input$stratum; treated <- input$treated
  denom.variance <- match.arg(denom.variance)
  solver <- match.arg(solver)
  .natural.flag(correct.sizes, "correct.sizes")
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(weight.by.size) || weight.by.size < 0 || weight.by.size > 1)
    stop("weight.by.size must lie between zero and one.")
  if (!scalar(runs) || runs < 1 || runs != floor(runs)) stop("runs must be a positive integer.")
  if (!is.null(seed) && (!scalar(seed) || seed < 0 || seed != floor(seed) || seed > .Machine$integer.max))
    stop("seed must be NULL or a nonnegative integer in the R seed range.")
  if (inherits(balance.formulas, "formula")) balance.formulas <- list(balance.formulas)
  if (!is.list(balance.formulas) || !length(balance.formulas) ||
      !all(vapply(balance.formulas, inherits, logical(1), "formula")))
    stop("balance.formulas must be a formula or a nonempty list of formulas.")
  if (!is.numeric(ratio) || !is.null(dim(ratio)) || !length(ratio) || any(!is.finite(ratio)) || any(ratio <= 0))
    stop("ratio must contain positive finite group ratios.")
  labels <- levels(z)
  if (length(ratio) == 1L) {
    targets <- stats::setNames(rep(as.numeric(ratio), length(labels)), labels)
    targets[treated] <- 1
  } else {
    if (is.null(names(ratio)) || anyNA(names(ratio)) || anyDuplicated(names(ratio)))
      stop("Multiple group ratios must be uniquely named by group.")
    if (setequal(names(ratio), labels[labels != treated])) ratio <- c(ratio, stats::setNames(1, treated))
    if (!setequal(names(ratio), labels)) stop("ratio names must identify all groups or all controls.")
    targets <- ratio[labels]
    if (targets[treated] != 1) stop("The treated group's target ratio must be one.")
  }
  requested <- round(outer(as.numeric(targets), as.numeric(input$counts[treated, ])))
  dimnames(requested) <- dimnames(input$counts)
  if (any(requested < 1)) stop("Every stratum must have positive target counts for each group.")
  if (any(requested > input$counts)) stop("The target ratio is infeasible in at least one stratum.")
  if ("stratum" %in% names(dat) && !identical(as.character(dat$stratum), as.character(stratum)))
    stop("dat$stratum disagrees with the supplied stratum vector.")
  design <- dat
  design$stratum <- stratum
  constraints <- .natural.backend(natstrat::generate_constraints(balance_formulas = balance.formulas,
    z = z, data = design, treated = treated, denom_variance = denom.variance,
    weight_by_size = weight.by.size))
  dots <- list(...)
  reserved <- c("z", "X", "st", "importances", "ratio", "treated", "seed", "runs", "solver", "correct_sizes", "q_s")
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) || any(names(dots) %in% reserved)))
    stop("Additional optimizer arguments must be named and cannot replace the primary design arguments.")
  backend <- do.call(natstrat::optimize_controls, c(list(z = z, X = constraints$X,
    st = stratum, importances = constraints$importances, ratio = unname(targets),
    treated = treated, seed = seed, runs = runs, solver = solver,
    correct_sizes = correct.sizes), dots))
  if (is.null(backend) || !is.logical(backend$selected) ||
      length(backend$selected) != nrow(dat) || anyNA(backend$selected))
    stop("The optimizer did not return a complete selection vector.")
  selected <- backend$selected
  counts <- table(z[selected], stratum[selected])
  if (any(counts < 1)) stop("A realized stratum lacks one group; revise the ratio or design.")
  ratios <- sweep(counts, 2, counts[treated, ], "/")
  proportional <- all(apply(ratios, 1, function(x) max(x) - min(x) < 1e-10))
  achieved <- all(abs(ratios - as.numeric(targets)) < 1e-10)
  list(data = dat[selected, , drop = FALSE], selected = selected,
       z = z[selected], stratum = stratum[selected], counts = counts,
       ratios = ratios, target.ratios = targets, target.counts = requested,
       proportional = proportional, target.ratio.achieved = achieved,
       treated = treated, constraints = constraints, backend = backend,
       seed = seed, solver = solver)
}
