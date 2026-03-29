#' @importFrom Formula as.Formula
#' @importFrom stats lm model.frame model.matrix model.response reformulate update
#' @importFrom np npcdistbw npcdist
#'
copula2sCOPEnp_fit <- function(F.formula, data, names.endo.reg, verbose){

  F.formula <- Formula::as.Formula(F.formula)

  #All endogenous regressors (continuous and discrete)
  endo.reg <- c(names.endo.reg$continuous, names.endo.reg$discrete) #discrete endo regressors are handled the same way
  #conditional CDF smooths the discrete distribution through X

  F.formula.main <- formula(F.formula, rhs = 1, lhs = 1)
  mf <- model.frame(F.formula.main, data = data)
  X.main <- model.matrix(F.formula.main, data = mf)

  endo.cols <- colnames(X.main)[colnames(X.main) %in% endo.reg]

  if(length(endo.cols) ==0)
    stop("No endogenous regressors found in the design matrix.", call. = FALSE)

  if(length(endo.col) < length(endo.reg))
    stop(paste0(
      "Bootstrap sample dropped at least one endogenous regressor. ",
      "This happened when a regressor becomes constant in the resampling."
    ), call. = FALSE)

  #exo column is everything that is non-intercept and non-endo col
  #used as conditioning variable X in the nonpara CDF

  exo.cols <- colnames(X.main)[!colnames(X.main) %in% c("Intercept", endo.cols)]

  if (length(exo.cols)==0)
    stop(paste0(
      "2sCOPEnp requires at least one exogenous regressor for the ",
      "nonparametric conditional CDF estimation.",
      "Please include exogenous control variables in the formula."
    ), call. = FALSE)

  np.data <- as.data.frame(X.main) #np functions need original data values not design matrix columns

  #first stage computing conditional CDF correction term nonparametrically (via helper) C_p
  #from equation 21 of Hu et al. 2025
  cop.term <- copula2sCOPEnp_correction(
    data = np.data,
    endo.cols = endo.cols,
    exo.cols = exo.cols,
    verbose = verbose
  )

  #second stage: augmented OLS. Adding the correction term to the structural model and estimate by OLS
  #using equation 20: Y = mu + sum_{k=1} ^ {K} ( P_{i,k} * alpha_k + beta' X_i + sum_{k=1}^{K} C_{i,pk} * gamma_k + epsilon_i

  f.main <- formula(mf)
  f.pcop <- reformulate(
    termlabels = c(".", colnames(cop.term)),
    response = NULL,
    intercept = TRUE
  )

  f.final <- update(old = f.main, new = f.pcop)

  return(lm(formula = f.final, data =cbind, cop.term))

}

