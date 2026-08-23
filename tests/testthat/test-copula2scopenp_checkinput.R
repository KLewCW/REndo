# tests only for the internal check helpers in f_copula2scopenp_checkinput.R
# these do 2sopeNP-specific input checks

# checkinput_copula2scopenp_formula_data --------------------------------------------
# test only what the aggregating method specifically checks for

test_that("formula_data fails for illegal specifications", {
  df <- fixture_copula_df()
  # Each case maps to a distinct rejected rule; regexp pins the specific branch.
  cases <- list(
    "continuous() special" = list(
      formula = y ~ x_num + x_fac | continuous(x_fac),
      data = df
    ),
    "discrete() special" = list(
      formula = y ~ x_num + x_fac | discrete(x_fac),
      data = df
    ),
    "interaction in rhs1" = list(formula = y ~ x_num * x_fac | x_fac, data = df),
    "no exogenous regr" = list(formula = y ~ x_num | x_num, data = df),
    "dot in rhs2 (all endo)" = list(formula = y ~ x_num + x_fac | ., data = df),
    "endo not in rhs1" = list(formula = y ~ x_num | x_fac, data = df)
  )
  check_errmsg(
    fn = checkinput_copula2scopenp_formula_data,
    cases = cases,
    regexp = c(
      "discrete\\(\\)/continuous\\(\\)",
      "discrete\\(\\)/continuous\\(\\)",
      "Interaction terms",
      "at least one exogenous",
      "at least one exogenous", # dot expands = everything endo = no exo
      "not in structural model"
    )
  )
})

test_that("formula_data works for valid formulas", {
  df <- fixture_copula_df()
  # to do log()
  df$x_pos <- abs(df$x_num) + 1

  cases <- list(
    # plain
    y ~ x_num + x_fac | x_num,
    # backtick exo
    y ~ x_num + `x space` | x_num,
    # backtick endo
    y ~ `x space` + x_num | `x space`,
    # transformation in response
    exp(y) ~ x_num + x_fac | x_num,
    # transformation in endo
    y ~ x_num + log(x_pos) | log(x_pos)
  )
  for (f in cases) {
    expect_null(checkinput_copula2scopenp_formula_data(formula = f, data = df))
  }

  # dot in rhs1
  # limit data to allowed types NA in data to test dot
  allowed.cols <- c("x_num", "x_fac", "x_ord", "y")
  expect_null(checkinput_copula2scopenp_formula_data(
    formula = y ~ . | x_num,
    data = df[, colnames(df) %in% allowed.cols]
  ))
})

# checkinput_copula2scopenp_npcdistbwargs -------------------------------------------
test_that("npcdistbwargs fails for invalid input", {
  cases <- list(
    "not a list" = list(npcdistbw.args = c(bwmethod = "cv.ls")),
    "unnamed" = list(npcdistbw.args = list("cv.ls")),
    "partially named" = list(npcdistbw.args = list(bwmethod = "cv.ls", "x")),
    "reserved xdat" = list(npcdistbw.args = list(xdat = 1)),
    "reserved ydat" = list(npcdistbw.args = list(ydat = 1))
  )
  check_errmsg(
    fn = checkinput_copula2scopenp_npcdistbwargs,
    cases = cases,
    regexp = c("plain list", "must be named", "must be named", "reserved", "reserved")
  )
})

test_that("npcdistbwargs warns for unknown arg name", {
  expect_warning(
    object = checkinput_copula2scopenp_npcdistbwargs(
      npcdistbw.args = list(abc = 1)
    ),
    regexp = "Unknown argument"
  )
})

test_that("npcdistbwargs works with defaults and legals args", {
  expect_null(checkinput_copula2scopenp_npcdistbwargs(npcdistbw.args = list()))
  expect_null(checkinput_copula2scopenp_npcdistbwargs(
    npcdistbw.args = list(bwmethod = "cv.ls", nmulti = 2)
  ))
})

# checkinput_copula2scopenp_bws -----------------------------------------------------
test_that("bws fails for illegal inputs", {
  bw <- structure(list(), class = "condbandwidth")
  cases <- list(
    "not a list" = list(bws = "x", labels.endo = "P"),
    "unnamed list" = list(bws = list(bw), labels.endo = "P"),
    "duplicate names" = list(bws = list(P = bw, P = bw), labels.endo = c("P", "Q")),
    "names not endo" = list(bws = list(Q = bw), labels.endo = "P"),
    "names not endo (multi)" = list(
      list(Q1 = bw, Q2 = bw),
      labels.endo = c("Q1", "P2")
    ),
    "wrong class" = list(bws = list(P = 1), labels.endo = "P")
  )
  check_errmsg(
    fn = checkinput_copula2scopenp_bws,
    cases = cases,
    regexp = c(
      "must be a list",
      "named list",
      "duplicate names",
      "named after",
      "named after",
      "condbandwidth"
    )
  )
})

test_that("bws works for NULL and valid list", {
  bw <- structure(list(), class = "condbandwidth")
  expect_null(checkinput_copula2scopenp_bws(bws = NULL, labels.endo = "P"))
  expect_null(checkinput_copula2scopenp_bws(
    bws = setNames(list(bw), "P"),
    labels.endo = "P"
  ))
})
