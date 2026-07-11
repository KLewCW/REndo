# Methods to apply expect_<X> to multiple test cases

# Expect error messages returned
check_errmsg <- function(fn, base.args = list(), cases, regexp) {
  if (length(regexp) == 1) {
    regexp <- rep(regexp, length(cases))
  }
  for (i in seq_along(cases)) {
    nm <- names(cases)[i]
    args <- c(base.args, cases[[i]])
    res <- do.call(what = fn, args = args)
    expect_true(length(res) > 0, label = nm)
    expect_match(paste(res, collapse = " "), regexp[i], info = nm)
  }
}

# Expect no error message returned
check_errmsg_empty <- function(fn, base.args = list(), cases) {
  for (nm in names(cases)) {
    expect_null(do.call(what = fn, args = c(base.args, cases[[nm]])))
  }
}

# check_clean_fit -------------------------------------------------------------------
check_clean_fit <- function(res) {
  expect_false(anyNA(fitted(res)))
  expect_false(anyNA(residuals(res)))
  expect_false(anyNA(coef(res)))
  expect_false(anyNA(vcov(res)))
  expect_false(anyNA(res$boots.params))
}
# Suppress boots warning -----------------
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

# check_param_recovery --------------------------
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

  cop.coefs  <- coef(res.lm)[aux.names]
  cop.matrix <- model.matrix(res.lm)[, aux.names, drop = FALSE]

  residuals.alt <- drop(residuals(res.lm) + cop.matrix %*% cop.coefs)
  fitted.alt    <- drop(fitted(res.lm) - cop.matrix %*% cop.coefs)

  expect_equal(residuals(res), residuals.alt)
  expect_equal(fitted.values(res), fitted.alt)
}
