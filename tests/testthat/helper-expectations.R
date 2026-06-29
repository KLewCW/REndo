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
