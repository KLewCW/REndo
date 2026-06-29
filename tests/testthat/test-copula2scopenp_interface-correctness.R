set.seed(42)

data("dataCopula2sCOPEnpCont")
data("dataCopula2sCOPEnpBi")
data("dataCopula2sCOPEnpMulti")


data.cont.pos <- local({
  df <- dataCopula2sCOPEnpCont
  df$P.pos <- abs(df$P) + 1
  df$X.pos <- abs(df$X) + 1
  return(df)
})

check_labelled_consistently <- function(res) {
  # aux terms in coefs
  expect_true(all(res$labels.pcop %in% names(coef(res))))
  # one aux term per endo
  expect_equal(length(res$labels.pcop), length(res$labels.endo))
  # one bw per endo, named after endo
  expect_setequal(names(res$bws), res$labels.endo)
}


# Invariance ---------------------------------------------------------------------------

test_that("Different sorted data produces same coefs", {
  # uses heuristic which should produce the exact same number
  res.sorted <- fit_2scopenp_fast(
    formula = y ~ P + X | P,
    data = data.cont.pos[order(data.cont.pos$y), ]
  )
  res.rev <- fit_2scopenp_fast(
    formula = y ~ P + X | P,
    data = data.cont.pos[rev(order(data.cont.pos$y)), ]
  )
  expect_equal(coef(res.sorted), coef(res.rev))
})

test_that("Order in formula does not impact coefs", {
  res.px <- fit_2scopenp_fast(formula = y ~ P + X | P, data = data.cont.pos)
  res.xp <- fit_2scopenp_fast(formula = y ~ X + P | P, data = data.cont.pos)
  expect_equal(sort(coef(res.px)), sort(coef(res.xp)))
})

test_that("Duplicate regressor has same result as single", {
  res.normal <- fit_2scopenp_fast(formula = y ~ P + X | P, data = data.cont.pos)
  res.dup <- fit_2scopenp_fast(formula = y ~ P + X + X | P, data = data.cont.pos)
  expect_equal(coef(res.dup), coef(res.normal))
})

test_that("Correctly respect intercept", {
  res.with <- fit_2scopenp_fast(formula = y ~ P + X | P, data = data.cont.pos)
  res.without <- fit_2scopenp_fast(formula = y ~ P + X - 1 | P, data = data.cont.pos)
  expect_true("(Intercept)" %in% names(coef(res.with)))
  expect_false("(Intercept)" %in% names(coef(res.without)))
  expect_equal(length(coef(res.with)), length(coef(res.without)) + 1)
})


# Internal consistency --------------------------------------------------------------

test_that("Formula edge cases: Label & residuals correct", {
  df.bt <- data.cont.pos
  df.bt$`P 1` <- df.bt$P

  df.multi.bt <- dataCopula2sCOPEnpMulti
  df.multi.bt$`P 1` <- df.multi.bt$P1

  cases <- list(
    "plain" = list(formula = y ~ P + X | P, data = data.cont.pos),
    "dot in rhs1" = list(formula = y ~ . | P, data = dataCopula2sCOPEnpCont),
    "backticks in endo" = list(formula = y ~ `P 1` + X | `P 1`, data = df.bt),
    "backticks in exo" = list(
      formula = y ~ P + `X 1` | P,
      data = transform(data.cont.pos, `X 1` = X, check.names = FALSE)
    ),
    "transformed endo" = list(
      formula = y ~ log(P.pos) + X | log(P.pos),
      data = data.cont.pos
    ),
    "multiple endo with backtick" = list(
      formula = y ~ `P 1` + P2 + X1 + X2 | `P 1` + P2,
      data = df.multi.bt
    )
  )

  for (nm in names(cases)) {
    res <- fit_2scopenp_fast(formula = cases[[nm]]$formula, data = cases[[nm]]$data)
    check_labelled_consistently(res)
    check_struct_residuals(res = res, aux.names = res$labels.pcop)
  }
})


test_that("Internals have correct shape", {
  res <- fit_2scopenp_fast(formula = y ~ P + X | P, data = data.cont.pos)
  n <- nrow(res$model)
  expect_length(fitted(res), n)
  expect_length(residuals(res), n)
  expect_equal(rownames(vcov(res)), names(coef(res)))
  expect_equal(colnames(vcov(res)), names(coef(res)))
})


# bws is used ----------------------------------------------------------------

test_that("Parameter bws is used as-is", {
  # When using heuristic, the result would be the same anyways. Therefore change it and
  # check that returns a different one
  res1 <- fit_2scopenp_fast(formula = y ~ P + X | P, data = data.cont.pos)

  bws.mod <- res1$bws
  for (nm in names(bws.mod)) {
    bws.mod[[nm]]$xbw <- bws.mod[[nm]]$xbw * 3
    bws.mod[[nm]]$ybw <- bws.mod[[nm]]$ybw * 3
  }

  res.perturbed <- fit_2scopenp_fast(
    formula = y ~ P + X | P,
    data = data.cont.pos,
    bws = bws.mod,
    npcdistbw.args = list(bandwidth.compute = FALSE)
  )

  # if bws were ignored the coefficients would be unchanged
  expect_false(isTRUE(all.equal(coef(res.perturbed), coef(res1))))
})


# Params recoveret ------------------------------------------------------------------

test_that("Recovery: Continuous endo (dataCopula2sCOPEnpCont)", {
  skip_on_cran()
  res <- fit_2scopenp_defaults(
    formula = y ~ P + X | P,
    data = dataCopula2sCOPEnpCont,
    npcdistbw.args = list(nmulti = 1)
  )
  check_param_recovery(res = res, true_vals = c("(Intercept)" = 1, P = 1, X = 2))
})

test_that("Recovery: Binary enod (dataCopula2sCOPEnpCont)", {
  skip_on_cran()
  res <- fit_2scopenp_defaults(
    formula = y ~ P + X | P,
    data = dataCopula2sCOPEnpBi,
    npcdistbw.args = list(nmulti = 1)
  )
  check_param_recovery(res = res, true_vals = c("(Intercept)" = 0, P = 1, X = 2))
})

test_that("Recovery: multiple endo (dataCopula2sCOPEnpMulti)", {
  skip_on_cran()
  res <- fit_2scopenp_defaults(
    formula = y ~ P1 + P2 + X1 + X2 | P1 + P2,
    data = dataCopula2sCOPEnpMulti,
    npcdistbw.args = list(nmulti = 1)
  )
  check_param_recovery(
    res = res,
    true_vals = c("(Intercept)" = 1, P1 = 1, P2 = 1, X1 = 2, X2 = -1)
  )
})
