skip_on_cran()
set.seed(42)

data("dataCopula2sCOPEnpCont")
data("dataCopula2sCOPEnpBi")
data("dataCopula2sCOPEnpMulti")


# Out of the box --------------------------------------------------------------------
test_that("works out of the box with default params", {
  res <- expect_message(copula2sCOPEnp(
    formula = y ~ P + X | P,
    data = dataCopula2sCOPEnpCont
  ))
  check_clean_fit(res)
})


# Formula -------------------------------------------------------------------------
test_that("works for edge-cases, incl print & print(summary())", {
  cases <- list(
    "dot in rhs1" = y ~ . | P,
    "transformed exo" = y ~ P + log(X_pos) | P,
    "transformed endo" = y ~ log(P_pos) + X | log(P_pos),
    "single endo" = y ~ P + X | P
  )

  df <- dataCopula2sCOPEnpCont
  df$X_pos <- abs(df$X) + 1
  df$P_pos <- abs(df$P) + 1

  for (nm in names(cases)) {
    res <- fit_2scopenp_fast(formula = cases[[nm]], data = df)
    check_clean_fit(res)
    # S3 methods
    expect_no_error(capture_output(print(res)))
    expect_no_error(capture_output(print(summary(res))))
  }
})

test_that("works with backticked variable names", {
  df <- dataCopula2sCOPEnpCont
  df$`P 1` <- df$P
  res <- fit_2scopenp_fast(formula = y ~ `P 1` + X | `P 1`, data = df)
  check_clean_fit(res)
})

test_that("works with multiple endo terms", {
  # only dataCopula2sCOPEnpMulti has >2 regressors
  # done use ordered factor columns because too slow
  res <- fit_2scopenp_fast(
    formula = y ~ X1 + X2 + P1 | P1 + X2,
    data = dataCopula2sCOPEnpMulti
  )
  check_clean_fit(res)
  expect_named(res$bws, c("P1", "X2"))
})


# Discrete endo ------------------------------------------------------------
test_that("works with a discrete endo", {
  res <- fit_2scopenp_fast(formula = y ~ P + X | P, data = dataCopula2sCOPEnpBi)
  check_clean_fit(res)
})

# Factor data -------------------------------------------------------------------------
test_that("works with ordered factors in in endo and exo", {
  res <- fit_2scopenp_fast(
    formula = y ~ P1 + P2 + X1 + X2 + X3 | P1 + P2,
    data = dataCopula2sCOPEnpMulti[1:250, ],
    # cannot run with heuristic
    npcdistbw.args = list(nmulti=1, tol=1, ftol=1)
  )
  check_clean_fit(res)
  expect_no_error(capture_output(print(res)))
  expect_no_error(capture_output(print(summary(res))))
})


# Tolerates irrelevant data ----------------------------------------------------------
test_that("works with NA in an irrelevant column", {
  df <- dataCopula2sCOPEnpCont
  df$abc <- rnorm(nrow(df))
  df$abc[1] <- NA
  res <- fit_2scopenp_fast(formula = y ~ P + X | P, data = df)
  check_clean_fit(res)
})


# Small sample ----------------------------------------------------------------------
test_that("works on a small sample", {
  small <- dataCopula2sCOPEnpCont[1:50, ]
  res <- fit_2scopenp_fast(formula = y ~ P + X | P, data = small)
  check_clean_fit(res)
})


# npcdistbw.args + bws: iterative refinement -----------------------------------------

test_that("works with a supplied bws and via iterative refinement", {
  res_loose <- suppress_lowboots_warning(
    copula2sCOPEnp(
      formula = y ~ P + X | P,
      data = dataCopula2sCOPEnpCont,
      npcdistbw.args = list(bwmethod = "normal-reference"),
      num.boots = 2,
      verbose = FALSE
    )
  )
  check_clean_fit(res_loose)

  # refine: pass previous bandwidths back as starting points
  res_refined <- suppress_lowboots_warning(
    copula2sCOPEnp(
      formula = y ~ P + X | P,
      data = dataCopula2sCOPEnpCont,
      bws = res_loose$bws,
      npcdistbw.args = list(nmulti = 1, ftol = 1, tol = 1),
      num.boots = 2,
      verbose = FALSE
    )
  )
  check_clean_fit(res_refined)
})


# verbose ---------------------------------------------------------------------------
test_that("verbose works", {
  expect_message(
    suppress_lowboots_warning(
      copula2sCOPEnp(
        formula = y ~ P + X | P,
        data = dataCopula2sCOPEnpCont,
        npcdistbw.args = list(bwmethod = "normal-reference"),
        num.boots = 2,
        verbose = TRUE
      )
    )
  )
  expect_no_message(
    fit_2scopenp_fast(formula = y ~ P + X | P, data = dataCopula2sCOPEnpCont)
  )
})
