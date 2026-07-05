# dont skip on cran: These are important tests! (and fast)
set.seed(42)

# checkinput_copulashared_data_basics ---------------------------------------
test_that("_data_basics rejects invalid data", {
  df <- fixture_copula_df()

  check_errmsg(
    checkinput_copulashared_data_basics,
    cases = list(
      "not a df" = list(data = list(1, 2)),
      "zero rows" = list(data = df[0, ]),
      "zero cols" = list(data = df[, 0])
    ),
    regexp = c(
      "must be a data.frame",
      "at least one row",
      "at least one column"
    )
  )

  # row + col messages
  expect_length(checkinput_copulashared_data_basics(data = df[0, 0]), 2)
})

test_that("data_basics accepts valid data", {
  # the data.frame itself is valid
  expect_null(checkinput_copulashared_data_basics(data = fixture_copula_df()))
})


# checkinput_copulashared_formula_basics --------------------------------------------
test_that("_formula_basics rejects invalid formulas", {
  check_errmsg(
    fn = checkinput_copulashared_formula_basics,
    cases = list(
      "not a formula" = list(formula = "y ~ x | z"),
      "one part" = list(formula = y ~ x),
      "more than 2" = list(formula = y ~ x | z | w),
      "multiple responses" = list(formula = y1 + y2 ~ x | z)
    ),
    regexp = c(
      "must be a formula",
      "two right",
      "two right",
      "exactly one response"
    )
  )
})

test_that("_formula_basics accepts valid formulas incl all kinds of edge cases", {
  # accepts normal formula, does not require Formula!
  cases <- list(
    "normal " = list(formula = y ~ x | z),
    "with dot" = list(formula = y ~ . | z),
    "backticks" = list(formula = `y space` ~ `x 1` | `z endo`),
    "transformations" = list(formula = log(y) ~ log(x) | exp(z)),
    "inline function" = list(formula = y ~ I(x^2) | z),
    "multi term" = list(formula = y ~ x1 + x2 | z1 + z2)
  )
  check_errmsg_empty(fn = checkinput_copulashared_formula_basics, cases = cases)
})
# checkinput_copulashared_vars_in_data -------------------------------------------
test_that("_vars_in_data reports missing variables", {
  check_errmsg(
    fn = checkinput_copulashared_vars_in_data,
    base.args = list(data = fixture_copula_df()),
    cases = list(
      "missing in rhs1" = list(F.formula = as.Formula(y ~ mis | x_num)),
      "missing in rhs2" = list(F.formula = as.Formula(y ~ x_num | mis)),
      "missing in response" = list(F.formula = as.Formula(mis ~ x_num | x_num)),
      "two missing listed" = list(F.formula = as.Formula(y ~ gone1 | gone2))
    ),
    regexp = c(
      "mis",
      "mis",
      "mis",
      "gone1.*gone2"
    )
  )
})
test_that("_vars_in_data works with valid, non-syntactic, dot", {
  df <- fixture_copula_df()
  check_errmsg_empty(
    fn = checkinput_copulashared_vars_in_data,
    base.args = list(data = df),
    cases = list(
      "non syntactic" = list(F.formula = as.Formula(y ~ `x space` | x_num)),
      "normal" = list(F.formula = as.Formula(y ~ x_num | x_fac)),
      "dot" = list(F.formula = as.Formula(y ~ . | x_num))
    )
  )
})

# checkinput_copulashared_response_not_in_rhs ------------------------------------------
test_that("_response_not_in_rhs fails if response in rhs1", {
  expect_match(
    checkinput_copulashared_response_not_in_rhs(
      F.formula = as.Formula(y ~ y + x_num),
      rhs1.terms = terms(
        as.Formula(y ~ y + x_num),
        lhs = 0,
        rhs = 1,
        data = fixture_copula_df()
      )
    ),
    regexp = "also appear"
  )

  # with transformation
  expect_match(
    checkinput_copulashared_response_not_in_rhs(
      F.formula = as.Formula(y ~ exp(y) + x_num),
      rhs1.terms = terms(
        as.Formula(y ~ log(y)),
        lhs = 0,
        rhs = 1,
        data = fixture_copula_df()
      )
    ),
    regexp = "also appear"
  )
})

test_that("_response_not_in_rhs works when response not in rhs", {
  F.formula <- as.Formula(y ~ x_num + x_fac | x_num)
  rhs.terms <- terms(F.formula, lhs = 0, rhs = 1, data = fixture_copula_df())
  expect_null(
    checkinput_copulashared_response_not_in_rhs(
      F.formula = F.formula,
      rhs1.terms = rhs.terms
    )
  )
})


# canonical_colname -------------------------------------------------------------------
test_that("canonical_colname normalizes labels to col names", {
  # plain
  expect_equal(canonical_colname(lab = "x_num"), "x_num")
  # backticks (space)
  expect_equal(canonical_colname(lab = "`x space`"), "x space")
  # cannot have space in label (requires backticks)
  # expect_equal(canonical_colname(lab = "x space"), "x space")
  # transformation
  expect_equal(canonical_colname(lab = "log(P)"), "log(P)")
})


# checkinput_copulashared_modelframe ---------------------------------------------------
test_that("_modelframe fails on NAs, illegal classes", {
  # NAs from transformation (also produces warning)
  expect_match(
    expect_warning(checkinput_copulashared_modelframe(
      F.formula = as.Formula(y ~ log(x_num - 1000)),
      data = fixture_copula_df(),
      allowed.classes = list("log(x_num - 1000)" = "numeric"),
      warn.low.card = FALSE
    )),
    regexp = "missing values"
  )

  check_errmsg(
    fn = checkinput_copulashared_modelframe,
    base.args = list(data = fixture_copula_df()),
    cases = list(
      "NA already in data" = list(
        F.formula = as.Formula(y ~ x_na),
        allowed.classes = list(x_na = "numeric"),
        warn.low.card = FALSE
      ),
      "given wrong class" = list(
        F.formula = as.Formula(y ~ x_fac),
        allowed.classes = list(x_fac = "numeric"),
        warn.low.card = FALSE
      )
    ),
    regexp = c(
      "missing values",
      "has class"
    )
  )
})

test_that("_modelframe warns for low-cardinality numerics", {
  df <- fixture_copula_df()
  expect_warning(
    checkinput_copulashared_modelframe(
      F.formula = as.Formula(y ~ x_lowcard),
      data = df,
      allowed.classes = list(x_lowcard = "numeric"),
      warn.low.card = TRUE
    ),
    regexp = "low cardinality"
  )
})

test_that("_modelframe works for valid data", {
  expect_null(
    checkinput_copulashared_modelframe(
      F.formula = as.Formula(y ~ x_num + x_fac + x_ord + `x space` + exp(x_num)),
      data = fixture_copula_df(),
      allowed.classes = list(
        x_num = "numeric",
        "exp(x_num)" = "numeric",
        "`x space`" = "numeric",
        x_fac = "factor",
        x_ord = "ordered"
      ),
      warn.low.card = TRUE
    )
  )
})
