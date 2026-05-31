#' @importFrom stats coef fitted residuals model.frame
new_rendo_copula2sCOPEnp <- function(
  call,
  F.formula,
  res.lm,
  boots.params,
  n.boots.attempted,
  n.boots.failed,
  names.endo.regs
) {
  return(.new_rendo_boots_degenerates_removed(
    # Stuff for rendo.boots.degenerates.removed
    call = call,
    F.formula = F.formula,
    mf = model.frame(res.lm),
    coefficients = coef(res.lm),
    names.main.coefs = row.names(boots.params),
    fitted.values = fitted(res.lm),
    residuals = resid(res.lm),
    boots.params = boots.params,
    n.boots.attempted = n.boots.attempted,
    n.boots.failed = n.boots.failed,

    # 2sCOPEnp-specific
    subclass = "rendo.copula.2sCOPE.np",
    names.endo.regs = names.endo.regs
  ))
}
