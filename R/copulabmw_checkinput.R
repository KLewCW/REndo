# Check formula (which needs data) and data
#
# Bare terms in RHS2: y ~ X + P | P
#
# - >= 1 endogenous regressor
# - >= 1 exogenous regressor
# - exogenous: any class
# - endogenous: numeric
# - wrongly specified continuous()/discrete() specials
#
#' @importFrom Formula as.Formula
#' @importFrom stats terms setNames
checkinput_copulabmw_formula_data <- function(formula, data) {
  err.msg <- checkinput_copulashared_formula_basics(formula = formula)
  if (length(err.msg) > 0) {
    return(err.msg)
  }

  F.formula <- as.Formula(formula)

  err.msg <- checkinput_copulashared_vars_in_data(F.formula = F.formula, data = data)
  if (length(err.msg) > 0) {
    return(err.msg)
  }

  # Expand dot against data
  # specials: check if user wrongly specified them
  specials <- c("discrete", "continuous")
  rhs1.terms <- terms(F.formula, lhs = 0, rhs = 1, data = data, specials = specials)
  rhs2.terms <- terms(F.formula, lhs = 0, rhs = 2, data = data, specials = specials)

  rhs1.labels <- labels(rhs1.terms)
  endo.labels <- labels(rhs2.terms)
  exo.labels <- rhs1.labels[!rhs1.labels %in% endo.labels]

  # Check after expanding dot
  err.msg <- checkinput_copulashared_response_not_in_rhs(
    F.formula = F.formula,
    rhs1.terms = rhs1.terms
  )
  if (length(err.msg) > 0) {
    return(err.msg)
  }

  # Catch if wrongly specified specials
  if (length(unlist(attr(rhs2.terms, "specials"))) > 0) {
    err.msg <- c(err.msg, "discrete()/continuous() are not supported for this method.")
  }

  # Need at least one of each: endogenous to model, exogenous for identification
  if (length(endo.labels) == 0) {
    err.msg <- c(err.msg, "At least one endogenous regressor (rhs2) is required.")
  }

  #equation 2.2 & assumption A4
  if (length(exo.labels) == 0) {
    err.msg <- c(
      err.msg,
      "At least one exogenous regressor is required for the first-stage regression of each endogenous regressor P on X."
    )
  }

  # Endo terms (RHS2) must also be in structural model (RHS1)
  not.in.rhs1 <- endo.labels[!endo.labels %in% rhs1.labels]
  if (length(not.in.rhs1) > 0) {
    # fmt: skip
    err.msg <- c(err.msg, paste0(
      "Endogenous regressor(s) not in structural model (rhs1): ",
      toString(not.in.rhs1), "."))
  }

  if (length(err.msg) > 0) {
    return(err.msg)
  }

  # endo: continuous
  # exo: anything - dont specify any classes for them
  allowed.classes <- setNames(rep(list("numeric"), length(endo.labels)), endo.labels)

  return(checkinput_copulashared_modelframe(
    F.formula = F.formula,
    data = data,
    allowed.classes = allowed.classes,
    labels.warn.low.card = endo.labels
  ))
}
