#' Adjusted Gaussian Copula Estimator (JAMS)
#'
#' @export
#' @importFrom stats coef model.response model.frame model.matrix
#'
copulaJAMS <- function(
  formula,
  data,
  cdf = c("adj.ecdf", "resc.ecdf", "ecdf", "kde"),
  num.boots = 1000,
  verbose = TRUE
) {
  cl <- match.call()

  # check_err_msg(checkinput_copulaJAMS_formula(formula))
  # check_err_msg(checkinput_copulaJAMS_data(data))
  # check_err_msg(checkinput_copulaJAMS_dataVSformula(data = data, formula = formula))
  # check_err_msg(checkinput_copulaJAMS_numboots(num.boots))
  # check_err_msg(checkinput_copulaJAMS_verbose(verbose))
  # check_err_msg(checkinput_copulaJAMS_cdf(cdf))

  cdf <- match.arg(cdf, choices = c("adj.ecdf", "resc.ecdf", "ecdf", "kde"))

  F.formula <- Formula::as.Formula(formula)
  f.main <- formula(F.formula, lhs = 1, rhs = 1)

  names.endo.regs <- formula_readout_special(
    F.formula = F.formula,
    name.special = "continuous",
    from.rhs = 2,
    params.as.chars.only = TRUE
  )

  #deriving exo regressor names from RHS1 - endo
  rhs.vars <- all.vars(formula(F.formula, rhs = 1, lhs = 0))
  names.exo.regs <- rhs.vars[!rhs.vars %in% names.endo.regs]

  #fitting the original data
  if (verbose) {
    message(
      "Fitting JAMS copula model for",
      length(names.endo.regs),
      "continuous endogenous regressor(s)."
    )
  }

  fit <- copulaJAMS_fit(
    f.main = f.main,
    data = data,
    names.endo.regs = names.endo.regs,
    names.exo.regs = names.exo.regs,
    cdf = cdf
  )

  # Bootstrapping ----------------------------------------------------------------------

  fn.fit.boots <- function(data.b) {
    return(
      copulaJAMS_fit(
        f.main = f.main,
        data = data.b,
        names.endo.regs = names.endo.regs,
        names.exo.regs = names.exo.regs,
        cdf = cdf
      )
    )
  }

  res.boots <- bootstrap_skip_degenerates(
    fn.fit = fn.fit.boots,
    data = data,
    num.boots = num.boots,
    coef.names = names(coef(fit)),
    verbose = verbose
  )

  # Structural residuals --------------------------------------------------------------

  l.fitted.resid <- copula_compute_structural_fitted_residuals(
    res.lm.aug = fit,
    names.aux.regs = grep("_cop$", names(coef(fit)), value = TRUE)
  )

  # Return object ----------------------------------------------------------------------

  return(new_rendo_copula2sCOPE(
    call = cl,
    F.formula = F.formula,
    res.lm.augmented = fit,
    fitted.values = l.fitted.resid$fitted.values,
    residuals = l.fitted.resid$residuals,
    boots.params = res.boots$boots.params,
    n.boots.attempted = res.boots$n.attempted,
    n.boots.failed = res.boots$n.failed,
    cdf = cdf,
    names.endo.regs = names.endo.regs
  ))
}
