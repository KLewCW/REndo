#'
#' @export
#' @importFrom stats coef
copulaBMW <- function(
    formula,
    data,
    cdf = c("ecdf", "adj.ecdf", "resc.ecdf", "kde"),
    num.boots = 1000,
    verbose = TRUE){

  cl <- match.call()

  check_err_msg(checkinput_copulaBMW_formula(formula))
  check_err_msg(checkinput_copulaBMW_data(data))
  check_err_msg(checkinput_copulaBMW_dataVSformula(data = data, formula = formula))
  check_err_msg(checkinput_copulaBMW_numboots(num.boots))
  check_err_msg(checkinput_copulaBMW_verbose(verbose))
  check_err_msg(checkinput_copulaBMW_cdf(cdf))

  cdf <- match.arg(cdf, choices = c("ecdf", "adj.ecdf", "resc.ecdf", "kde"))

  F.formula <- Formula::as.Formula(formula)
  names.endo.regs <- formula_readout_special(
    F.formula = F.formula,
    name.special = "continuous",
    from.rhs = 2,
    params.as.chars.only = TRUE
  )

  rhs1.vars <- all.vars(formula(F.formula, rhs = 1, lhs = 0))
  exo.vars <- rhs1.vars[!rhs1.vars %in% names.endo.regs]

  if(length(exo.vars) == 0){
    stop("No exogenous regressors were found. BMW method requires at least one",
         "exogeous regressor for the first-stage regression of each endogenous regressor ",
         "P on X.", call. = FALSE) #equation 2.2 & assumption A4
  }

  if(verbose){
    message("Fitting BMW copula model with",
            length(names.endo.regs),
            "endogenous regressors.")
  }

  fit <- copulaBMW_fit(
    F.formula = F.formula,
    data = data,
    names.endo.regs = names.endo.regs,
    cdf = cdf
  )

  fn.fit.boots <- function(data.b){
    return(copulaBMW_fit(F.formula = F.formula, data = data.b, cdf = cdf))
  }

  res.boots <- bootstrap_skip_degenerates(
    fn.fit     = fn.fit.boots,
    data       = data,
    num.boots  = num.boots,
    coef.names = names(coef(fit)),
    verbose    = verbose
  )

  return(new_rendo_copulaBMW(
    call              = cl,
    F.formula         = F.formula,
    res.lm            = fit,
    boots.params      = res.boots$boots.params,
    n.boots.attempted = res.boots$n.attempted,
    n.boots.failed    = res.boots$n.failed,
    cdf               = cdf,
    names.endo.regs   = names.endo.regs
  ))
}
