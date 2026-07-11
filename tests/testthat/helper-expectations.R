# check_errmsg() + check_errmsg_empty -------------------------------------------------
# Mostly for testing the internal methods that return error messages (char vecs) and
# do not stop()

# Methods to apply expect_<X> to multiple test cases
# Expect error messages returned
check_errmsg <- function(fn, base.args = list(), cases, regexp) {
  if (length(regexp) == 1) {
    regexp <- rep(regexp, length(cases))
  }
  # for printing labels on failure
  names(regexp) <- names(cases)
  case.i <- 0

  run_cases(
    fn = fn,
    base.args = base.args,
    cases = cases,
    check = function(res) {
      # `<<-` to manipulate up the call stack (counter in enclosing method)
      case.i <<- case.i + 1
      expect_true(length(res) > 0)
      expect_match(paste(res, collapse = " "), regexp[[case.i]])
    }
  )
}

# Expect no error message returned
check_errmsg_empty <- function(fn, base.args = list(), cases) {
  run_cases(
    fn = fn,
    base.args = base.args,
    cases = cases,
    check = expect_null
  )
}

# check_clean_fit -------------------------------------------------------------------
# No NAs anywhere. Basic sanity checks. Mostly used for smoke tests.
check_clean_fit <- function(res) {
  expect_false(anyNA(fitted(res)))
  expect_false(anyNA(residuals(res)))
  expect_false(anyNA(coef(res)))
  expect_false(anyNA(vcov(res)))
  expect_false(anyNA(res$boots.params))
}

# run_cases ------------------------------------------------------------------
# run `fn` method with `cases` + `base.args` and then apply `check` as a test
#' @importFrom utils modifyList
run_cases <- function(fn, base.args = list(), cases, check) {
  results <- lapply(names(cases), function(nm) {
    res <- do.call(what = fn, args = modifyList(base.args, cases[[nm]]))
    withCallingHandlers(
      check(res),
      # catch thrown testthat expectation failures
      expectation_failure = function(e) {
        e$message <- paste0("[case: ", nm, "] ", e$message)
        stop(e)
      }
    )
    return(res)
  })
  names(results) <- names(cases)
  return(invisible(results))
}

# Suppress boots warning -----------------------------------------------------------
# Suppress 1000 num.boots warning but lets all other warnings propagate
suppress_lowboots_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(pattern = "recommended to run 1000", x = conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# check_param_recovery -------------------------------------------------------------
# parameter recovery test
check_param_recovery <- function(res, true_vals) {
  coefs <- coef(res)
  ses <- sqrt(diag(vcov(res)))
  for (nm in names(true_vals)) {
    diff <- abs(coefs[nm] - true_vals[nm])
    expect_true(
      object = diff < 2 * ses[nm],
      info = sprintf(
        "%s: est=%.3f, true=%.3f, 2*SE=%.3f",
        nm,
        coefs[nm],
        true_vals[nm],
        2 * ses[nm]
      )
    )
  }
}

# check_struct_residuals ----------------------------------------------------------
#' @importFrom stats residuals fitted coef
check_struct_residuals <- function(res, aux.names) {
  # Alternative route: Remove cop contribution from augmented fit
  res.lm <- res$res.lm.augmented

  cop.coefs <- coef(res.lm)[aux.names]
  cop.matrix <- model.matrix(res.lm)[, aux.names, drop = FALSE]

  residuals.alt <- drop(residuals(res.lm) + cop.matrix %*% cop.coefs)
  fitted.alt <- drop(fitted(res.lm) - cop.matrix %*% cop.coefs)

  expect_equal(residuals(res), residuals.alt)
  expect_equal(fitted.values(res), fitted.alt)
}
