# Structure: y ~ X + P | P
# - >= 1 endogenous regressors
# - >= 1 exogenous regressors
# - all regressors and response have to be numeric
# - Interactions: Block (because uncertain whether really allowed, eg interaction with endo var)
#
#' @importFrom Formula as.Formula
#' @importFrom stats terms setNames formula
checkinput_copulabayes_formula_data <- function(formula, data) {
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

  # Wrongly specified specials
  if (length(unlist(attr(rhs2.terms, "specials"))) > 0) {
    err.msg <- c(err.msg, "discrete()/continuous() are not supported for this method.")
  }

  # Need at least one exo & endo
  if (length(endo.labels) == 0) {
    err.msg <- c(err.msg, "At least one endogenous regressor (rhs2) is required.")
  }
  if (length(exo.labels) == 0) {
    err.msg <- c(err.msg, "At least one exogenous regressor is required.")
  }

  # Endogenous terms (RHS2) must also be in structural model (RHS1)
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

  # No interactions
  if (any(attr(rhs1.terms, "order") > 1)) {
    return("Interaction terms are not supported for this method.")
  }

  # endo and exo: all numeric
  allowed.classes <- setNames(rep(list("numeric"), length(rhs1.labels)), rhs1.labels)

  return(checkinput_copulashared_modelframe(
    F.formula = F.formula,
    data = data,
    allowed.classes = allowed.classes,
    labels.warn.low.card = endo.labels
  ))
}


checkinput_copulabayes_numiterations_burnin_thin <- function(
  num.iterations,
  burnin,
  thin
) {
  err.msg <- c(
    checkinputhelper_singlepositivewholenumeric(
      num.param = num.iterations,
      parameter.name = "num.iterations",
      min.num = 1
    ),
    checkinputhelper_singlepositivewholenumeric(
      num.param = burnin,
      parameter.name = "burnin",
      min.num = 0
    ),
    checkinputhelper_singlepositivewholenumeric(
      num.param = thin,
      parameter.name = "thin",
      min.num = 1
    )
  )
  if (length(err.msg) > 0) {
    return(err.msg)
  }

  # cannot burn in longer than full iterations
  if (burnin >= num.iterations) {
    return("Parameter 'burnin' must be smaller than 'num.iterations'.")
  }

  num.remaining <- floor((num.iterations - burnin) / thin)

  # Need at least 2 observations for operations like sd()
  if (num.remaining < 2) {
    # fmt: skip
    return(paste0("Only ", num.remaining, " draw(s) will be retained after thinning. Need at least 2."))
  }

  if (num.remaining < 1000) {
    # fmt: skip
    warning(
      "It is recommended to retain 1000 or more draws after burnin and thinning, ",
      "but the given parameters retain ", num.remaining, ".",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  return(c())
}
