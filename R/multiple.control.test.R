multiple.control.test <- function(y = NULL, group = NULL, means = NULL,
                                  se = NULL, n = NULL, treated = NULL,
                                  controls = NULL, alpha = 0.05,
                                  alternative = c("greater", "less"),
                                  na.rm = FALSE) {
  alternative <- match.arg(alternative)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1 || 1 - alpha >= 1 || 1 - alpha <= 0)
    stop("alpha and its confidence level 1-alpha must be strictly between zero and one.")
  groups <- .contrast.groups(y, group, means, se, n, na.rm)
  labels <- groups$labels
  if (length(labels) != 3L) stop("Exactly one treated group and two control groups are required.")
  if (is.null(treated)) {
    if (groups$raw) stop("Specify the treated group label for raw data.")
    treated <- labels[1]
  }
  if (length(treated) != 1L || is.na(treated) || !as.character(treated) %in% labels)
    stop("treated must identify one of the three group labels.")
  treated <- as.character(treated)
  if (is.null(controls)) controls <- labels[labels != treated]
  controls <- as.character(controls)
  if (length(controls) != 2L || anyNA(controls) || anyDuplicated(controls) ||
      !setequal(controls, labels[labels != treated]))
    stop("controls must name the two distinct remaining groups, in the desired order.")
  order <- c(treated, controls)
  coefficients <- matrix(c(1, -0.5, -0.5,
                           1, -1, 0,
                           1, 0, -1,
                           1, -2, 1,
                           1, 1, -2), nrow = 5, byrow = TRUE,
                         dimnames = list(paste0("H", 1:5), order))
  tests <- data.frame(hypothesis = rownames(coefficients),
                      stage = c(1L, 2L, 2L, 3L, 3L),
                      null = if (alternative == "greater") "contrast <= 0" else "contrast >= 0",
                      estimate = NA_real_, SE = NA_real_, df = NA_real_,
                      statistic = NA_real_, p.value = NA_real_, log.p.value = NA_real_,
                      tested = FALSE, rejected = NA, status = "skipped",
                      row.names = NULL)
  components <- list()
  run <- function(i) {
    fit <- contrast.test(contrast = coefficients[i, ], means = groups$means,
                         se = groups$se, n = groups$n, alternative = alternative,
                         conf.level = 1 - alpha)
    components[[rownames(coefficients)[i]]] <<- fit
    tests$estimate[i] <<- unname(fit$estimate)
    tests$SE[i] <<- fit$stderr
    tests$df[i] <<- unname(fit$parameter)
    tests$statistic[i] <<- unname(fit$statistic)
    tests$p.value[i] <<- fit$p.value
    tests$log.p.value[i] <<- fit$log.p.value
    tests$tested[i] <<- TRUE
    tests$rejected[i] <<- fit$p.value <= alpha
    tests$status[i] <<- if (fit$p.value <= alpha) "rejected" else "not rejected"
  }
  run(1L)
  if (tests$rejected[1]) {
    run(2L)
    run(3L)
    if (all(tests$rejected[2:3])) {
      run(4L)
      run(5L)
    }
  }
  rejected <- function(i) isTRUE(tests$rejected[i])
  conclusions <- c(treated.vs.average = rejected(1L),
                   treated.vs.control1 = rejected(2L),
                   treated.vs.control2 = rejected(3L),
                   treated.vs.both = rejected(2L) && rejected(3L),
                   separation.exceeds.control.difference = rejected(4L) && rejected(5L))
  list(tests = tests, coefficients = coefficients,
       groups = stats::setNames(order, c("treated", "control1", "control2")),
       group.summary = groups$summary, conclusions = conclusions,
       stage.reached = max(tests$stage[tests$tested]),
       component.tests = components, alpha = alpha, alternative = alternative,
       method = "Ordered one-sided Welch contrasts with two control groups",
       fwer.control = "Approximate: the ordering uses Welch component tests")
}
