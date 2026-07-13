# For each method:
#   fit_<method>_fast
#   fit_<method>_cv

# copula2sCOPEnp -----------------------------------------------------------------------

fit_2scopenp_fast <- function(
  formula,
  data,
  num.boots = 2,
  verbose = FALSE,
  npcdistbw.args = list(bwmethod = "normal-reference"),
  bws = NULL
) {
  return(suppress_lowboots_warning(
    copula2sCOPEnp(
      formula = formula,
      data = data,
      num.boots = num.boots,
      verbose = verbose,
      npcdistbw.args = npcdistbw.args,
      bws = bws
    )
  ))
}

fit_2scopenp_defaults <- function(
  formula,
  data,
  num.boots = 1000,
  verbose = FALSE,
  npcdistbw.args = list(),
  bws = NULL
) {
  return(suppress_lowboots_warning(
    copula2sCOPEnp(
      formula = formula,
      data = data,
      num.boots = num.boots,
      verbose = verbose,
      npcdistbw.args = npcdistbw.args,
      bws = bws
    )
  ))
}
