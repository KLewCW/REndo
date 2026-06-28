# Fixtures used in tests

# copula methods  ----------------------------------------------------------------------
fixture_copula_df <- function(n = 100) {
  return(data.frame(
    y = rnorm(n),
    x_num = rnorm(n),
    `x space` = rnorm(n), # non-syntactic name
    x_string = rep(c("a", "b"), length.out = n),
    x_fac = factor(rep(c("a", "b"), length.out = n)),
    x_ord = ordered(rep(c("lo", "hi"), length.out = n), levels = c("lo", "hi")),
    x_lowcard = rep(1:2, length.out = n),
    x_na = c(NA, rnorm(n - 1)),
    check.names = FALSE, # dont need to be syntactically valid
    stringsAsFactors = FALSE
  ))
}

