#'
#' @importFrom np npcdistbw npcdist
#' @importFrom stats qnorm

#computing the nonparametric conditional CDF correction terms

#For each endogenous regressor P_k, we have to estimate the conditional CDF
#F hat_(P_k |X) nonparametrically using Kernel methods (Hu et al 2025, mentions mostly Nadaraya-Watson (NW) kernel
#before eq. 10) Stage 1 from table 3
#Then the normal quantile transformation is applied to get the copula correction term
#' @importFrom np npcdistbw npcdist
copula2sCOPEnp_bandwidth <- function(y.data, x.data){
  h <- np::npcdistbw(
    xdat = x.data,
    ydat = y.data
  )
  return(h)
}
copula2sCOPEnp_correction <- function(data, endo.cols, exo.cols, verbose, bws = NULL) {
  n <- nrow(data)
  res <- matrix(NA_real_, nrow = n, ncol = length(endo.cols))
  colnames(res) <- paste0(endo.cols, "_cop")

  for (k in seq_along(endo.cols)) {
    p.var <- endo.cols[k]

    if (verbose) {
      message(
        "Computing conditional CDF for endogenous regressor ' ",
        p.var,
        " ' (",
        k,
        " of  ",
        length(endo.cols),
        ") ..."
      )
    }

    #extracting the endogenous regressor as a data.frame for the np functions

    y.data <- data[, p.var, drop = FALSE]

    #extracting the exogenous regressors as data.frame
    #npcdistbw: supposed to handle both cts and discrete (or combination of both) X automatically through a
    #generalised product kernel (Li and Racine 2008)

    x.data <- data[, exo.cols, drop = FALSE]

    #selecting bandwidth (h) for conditional CDF F hat_(P_k |X) through cross validation

    if (is.null(bws)){
      h <- copula2sCOPEnp_bandwidth(
        y.data = y.data,
        x.data = x.data
      )
    } else{
      h <- bws[[k]] # reusing pre-computed bandwidths
    }

    #nonpara conditional CDF estimate

    cdf.fit <- np::npcdist(bws = h, newdata = data.frame(y.data, x.data))

    #calculating F hat_ (P_{i,k} | X_i) for each observation

    conditional.cdf <- cdf.fit$condist

    #applying normal quantle transformation:  C_{i,pk} = phi^{-1} (F hat_ (P_k | X))
    res[, k] <- qnorm(conditional.cdf) #table 3 stage 1 and equation 21 from Hu et al. 2025
  }

  return(res)
}
