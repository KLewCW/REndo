skip_on_cran()
set.seed(42)

data("dataCopBMW")
data("dataCopBMWMultiEndo")

all.cdfs <- c("ecdf", "adj.ecdf", "resc.ecdf", "kde")


# Out of the box --------------------------------------------------------------------
test_that("Works out of the box", {
  res <- expect_message(copulaBMW(
    formula = y ~ X + P | P,
    data = dataCopBMW
  ))
  check_clean_fit(res)
})


# specs x cdf ----------------------------------------------------------------------
# Specs vary along:
# - intercept: w/ and w/o
# - endo: 1 vs multiple
#
# All combinations: cross each spec with all 4 cdf options
test_that("Works for all specs and cdf combinations", {
  specs <- list(
    "with intercept, single endo" = list(
      formula = y ~ X + P | P,
      data = dataCopBMW
    ),
    "no intercept, single endo" = list(
      formula = y ~ X + P - 1 | P,
      data = dataCopBMW
    ),
    "with intercept, multiple endo" = list(
      formula = y ~ X + P1 + P2 | P1 + P2,
      data = dataCopBMWMultiEndo
    ),
    "no intercept, multiple endo" = list(
      formula = y ~ X + P1 + P2 - 1 | P1 + P2,
      data = dataCopBMWMultiEndo
    )
  )

  for (cdf in all.cdfs) {
    run_cases(
      fn = fit_bmw_fast,
      base.args = list(cdf = cdf),
      cases = specs,
      check = check_clean_fit
    )
  }
})


# Formula -------------------------------------------------------------------------
test_that("Works for formula edge-cases, incl basic s3 (print, summary)", {
  df <- dataCopBMW
  df$X_pos <- abs(df$X) + 1
  df$P_pos <- abs(df$P) + 1
  df$x_fac <- factor(rep(c("a", "b", "c"), length.out = nrow(df)))
  df$x_bin <- rep(0:1, length.out = nrow(df))
  df$`P 1` <- df$P
  df$`X 1` <- df$X

  cases <- list(
    "dot in rhs1" = list(formula = y ~ . | P, data = dataCopBMW),
    "backticks in formula" = list(formula = y ~ `P 1` + `X 1` | `P 1`),
    "transformed exo" = list(formula = y ~ log(X_pos) + P | P),
    "transformed endo" = list(formula = y ~ log(P_pos) + X | log(P_pos)),
    "factor exo as only" = list(formula = y ~ x_fac + P | P),
    "factor exo with cont" = list(formula = y ~ X + x_fac + P | P),
    # other form of discrete - but leads to low cardinality warning
    "binary exo" = list(formula = y ~ x_bin + P | P)
  )

  results <- run_cases(
    fn = fit_bmw_fast,
    base.args = list(data = df),
    cases = cases,
    check = check_clean_fit
  )

  # S3 methods on each edge case
  for (nm in names(results)) {
    expect_no_error(capture_output(print(results[[nm]])))
    expect_no_error(capture_output(print(summary(results[[nm]]))))
  }
})


# Tolerate irrelevant data ----------------------------------------------------------
test_that("Works with NA in irrelevant col", {
  df <- dataCopBMW
  df$abc <- rnorm(nrow(df))
  df$abc[1] <- NA
  res <- fit_bmw_fast(formula = y ~ X + P | P, data = df)
  check_clean_fit(res)
})


# Small sample ----------------------------------------------------------------------
test_that("Works on small sample", {
  small <- dataCopBMW[1:50, ]
  res <- fit_bmw_fast(formula = y ~ X + P | P, data = small)
  check_clean_fit(res)
})
