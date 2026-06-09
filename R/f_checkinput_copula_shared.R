checkinput_copulashared_cdf <- function(cdf, allowed.cdf) {
  return(checkinputhelper_choice(
    choice = cdf,
    allowed = allowed.cdf,
    param.name = "cdf"
  ))
}

checkinput_copulashared_numboots <- function(num.boots) {
  return(checkinputhelper_numboots(num.boots))
}

checkinput_copulashared_data <- function(data) {
  return(.checkinputhelper_data_basicstructure(data))
}

checkinput_copulashared_verbose <- function(verbose) {
  checkinputhelper_single_logical(logical = verbose, param.name = "verbose")
}


#' @importFrom Formula as.Formula
checkinput_copulashared_dataVSformula <- function(data, formula) {
  F.formula <- Formula::as.Formula(formula)

  names.cols.endo <- all.vars(terms(F.formula, rhs = 2, lhs = 0))

  err.msg <- .checkinputhelper_dataVSformula_basicstructure(
    formula = F.formula,
    data = data,
    rhs.rel.regr = c(1, 2),
    num.only.cols = names.cols.endo
  )

  # Still relevant??
  err.msg <- c(
    err.msg,
    checkinputhelper_data_notnamed(
      formula = F.formula,
      data = data,
      forbidden.colname = "Pstar"
    )
  )

  return(err.msg)
}


checkinput_copulashared_formula <- function(formula) {
  err.msg <- .checkinputhelper_formula_basicstructure(formula = formula)
  if (length(err.msg) > 0) {
    return(err.msg)
  }

  F.formula <- Formula::as.Formula(formula)

  # Check if exactly 2 part formula ---------------------------------------------------

  # check to see if the formula inputed has 2 RHS
  if (length(F.formula)[2] < 2) {
    err.msg <- c(
      err.msg,
      "Please specify endogenous regressors on a second right-hand side using '|' (e.g. y ~ X + P | continuous(P))."
    )
  }

  if (length(F.formula)[2] > 2) {
    err.msg <- c(
      err.msg,
      "Please specify only a two-part formula using only one separator '|'."
    )
  }

  if (length(err.msg) > 0) {
    return(err.msg)
  }

  # Check if endo are in structural model -----------------------------------------------

  #for the variables
  rhs1.vars <- all.vars(formula(F.formula, rhs = 1, lhs = 0))
  rhs2.vars <- all.vars(formula(F.formula, rhs = 2, lhs = 0))

  #RHS2 vars must also be in RHS1
  if (!all(rhs2.vars %in% rhs1.vars)) {
    err.msg <- c(
      err.msg,
      "Please specify every endogenous regressor also in the first right-hand side of the formula."
    )
  }

  # Checks related to RHS2 continuous() ------------------------------------------------

  # Checks for specials function `continuous()`
  names.vars.continuous <- formula_readout_special(
    F.formula = F.formula,
    name.special = "continuous",
    from.rhs = 2,
    params.as.chars.only = TRUE
  )

  # At least 1
  if (length(names.vars.continuous) == 0) {
    err.msg <- c(
      err.msg,
      "The method requires at least one continuous endogenous regressor, specified using `continuous()`."
    )
  }

  # Every RHS2 variable must be wrapped in continuous()
  rhs2.calls <- as.list(attr(terms(F.formula, rhs = 2, lhs = 0), "variables"))[-1]
  if (
    !all(vapply(
      rhs2.calls,
      function(e) {
        is.call(e) && identical(e[[1]], as.name("continuous"))
      },
      logical(1)
    ))
  ) {
    err.msg <- c(
      err.msg,
      "Please wrap every endogenous regressor on the second RHS in continuous()."
    )
  }

  # Variable names in continuous() must match exactly as they appear in RHS1
  # (including any transformations such as log(x) or I(x^2))
  rhs1.labels <- labels(terms(F.formula, rhs = 1, lhs = 0))
  if (!all(names.vars.continuous %in% rhs1.labels)) {
    err.msg <- c(
      err.msg,
      "Please name every endogenous regressor exactly as in the main model, including transformations if any."
    )
  }

  # Check if continuous() in LHS or first RHS ------------------------------------------

  num.specials.rhs1 <- sum(sapply(
    attr(terms(F.formula, rhs = 1, lhs = 0, specials = "continuous"), "specials"),
    length
  ))

  num.specials.lhs <- sum(sapply(
    attr(terms(F.formula, rhs = 0, lhs = 1, specials = "continuous"), "specials"),
    length
  ))

  if (num.specials.rhs1 > 0) {
    err.msg <- c(err.msg, "No endogenous regressors should be on the first RHS.")
  }

  if (num.specials.lhs > 0) {
    err.msg <- c(err.msg, "No endogenous regressors should be on the LHS.")
  }

  return(err.msg)
}


# NEW INPUT CHECKS ------------------------------------------------------------------


# Basic structure of data
# - is data.frame
# - has row & cols (>0)
#
checkinput_copulashared_data_basics <- function(data) {
  err.msg <- c()

  if (!is.data.frame(data)) {
    return("Parameter 'data' must be a data.frame.")
  }

  if (nrow(data) == 0) {
    err.msg <- c(err.msg, "Parameter 'data' must have at least one row.")
  }

  if (ncol(data) == 0) {
    err.msg <- c(err.msg, "Parameter 'data' must have at least one column.")
  }

  return(err.msg)
}


# Basic structure of the formula
#  - formula
#  - two part
#  - single response
#
#' @importFrom Formula as.Formula
checkinput_copulashared_formula_basics <- function(formula) {
  if (!inherits(formula, "formula")) {
    return("Parameter 'formula' must be a formula.")
  }

  # Report coercion error
  F.formula <- tryCatch(as.Formula(formula), error = function(e) conditionMessage(e))
  if (!inherits(F.formula, "Formula")) {
    return(paste0(
      "Parameter 'formula' could not be interpreted as a Formula: ",
      F.formula
    ))
  }

  parts <- length(F.formula)
  if (parts[1] != 1 || parts[2] != 2) {
    return(
      "Parameter 'formula' must have one left-hand side and two right-hand parts (lhs ~ rhs1 | rhs2)."
    )
  }

  # Check multiple responses (y1 + y2 ~ )
  if (length(all.vars(formula(F.formula, lhs = 1, rhs = 0))) != 1) {
    return("Parameter 'formula' must have exactly one response variable.")
  }

  return(c())
}

# All formula variables exist in the data
checkinput_copulashared_vars_in_data <- function(F.formula, data) {
  # Use unexpanded formula and not terms object: Check what the user actually specified

  # Also correct if there are non-syntactic names
  # all.vars() returns pure strings without backticks, like they are in colnames
  # (eg all.vars(~ `X space`) = "X space" = colnames("X space"))
  missing.vars <- setdiff(all.vars(F.formula), colnames(data))
  if (length(missing.vars) > 0) {
    return(paste0(
      "Variable(s) in parameter 'formula' not found in parameter 'data': ",
      paste0(missing.vars, collapse = ", "),
      "."
    ))
  }

  return(c())
}

# LHS is not also in RHS1
checkinput_copulashared_response_not_in_rhs <- function(F.formula, rhs1.terms) {
  # needs terms object of expanded terms object
  # use all.vars() and not labels as its about the variable ( y ~ log(y) is illegal)
  lhs.vars <- all.vars(formula(F.formula, lhs = 1, rhs = 0))
  rhs.vars <- all.vars(rhs1.terms)

  lhs.in.rhs <- intersect(lhs.vars, rhs.vars)
  if (length(lhs.in.rhs) > 0) {
    # fmt: skip
    return(paste0(
      "Response variable(s) in parameter 'formula' also appear on the right-hand side: ",
      toString(lhs.in.rhs),"."
    ))
  }

  return(c())
}

# Normalizes a term label to column names (as in model.frame())
canonical_colname <- function(lab) {
  # strips backticks from non-syntactic names
  return(deparse1(parse(text = lab)[[1]]))
}

# Build model.frame and verify column classes
# - build from lhs ~ rhs1
# - expose failed builds to user
# - na.fail to catch generated NAs
# - check for each (formula) label / term, if the resulting regressor (colum) is
#    of an allowed class
# - warn about low cardinality (<=10) numeric variables
#
#' @importFrom stats model.frame na.fail .MFclass
checkinput_copulashared_modelframe <- function(F.formula, data, allowed.classes) {
  # Expose build failures to user (bad transformations, generated NAs,...)
  mf <- tryCatch(
    # Making model.frame with rhs=1 also verifies that rhs=2 regressors are fine
    model.frame(formula(F.formula, lhs = 1, rhs = 1), data = data, na.action = na.fail),
    error = function(e) conditionMessage(e)
  )
  if (!is.data.frame(mf)) {
    return(paste0(
      "The model frame (response ~ structural regressors) could not be built: ",
      mf
    ))
  }

  err.msg <- c()
  mf.cols <- colnames(mf)

  # For each column: Check if is allowed class
  for (lab in names(allowed.classes)) {
    col <- canonical_colname(lab)

    # Should not happen: Previous check should verify already that all formula labels
    # are in data
    stopifnot(col %in% mf.cols)

    cls <- .MFclass(mf[[col]])

    if (!cls %in% allowed.classes[[lab]]) {
      # fmt: skip
      err.msg <- c(err.msg, paste0(
        "Regressor '", lab, "' has class '", cls,
        "' but must be one of: ", toString(allowed.classes[[lab]]), "."))
    }
  }

  # warn about low-cardinality numeric variables
  is.low.card <- sapply(mf, function(x){
    # no NAs at this point
    is.numeric(x) && length(unique(x)) <= 10
    })

  if(any(is.low.card)){

    low.card.vars <- names(which(is.low.card))

    warning(
      "The following numeric regressors have low cardinality (<= 10 distinct values) and may be non-continuous: ",
      toString(low.card.vars),
      "\nTreating them as numeric instead of (ordered) factors may yield wrong results as the method depend on correctly classed data. ",
      "Consider converting them with `factor()` or `ordered()`, if appropriate.",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  return(err.msg)
}


