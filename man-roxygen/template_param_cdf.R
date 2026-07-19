#' @param cdf Character string specifying the method used to estimate the
#'   marginal distribution functions of the endogenous regressors.
#'   One of the following:
#' \describe{
#'     \item{\code{"adj.ecdf"}}{Adjusted empirical CDF with midrank correction
#'       (Liengaard et al., 2024). Keeps values strictly inside (0,1).}
#'     \item{\code{"resc.ecdf"}}{Rescaled empirical CDF via
#'       \code{copula::pobs} (Qian et al., 2024).}
#'      \item{\code{"ecdf"}}{<%= if (exists("ecdf.text")) ecdf.text else "Empirical CDF with boundary replacement (Becker et al., 2022)." %>}
#'     \item{\code{"kde"}}{Integral of a density estimator via \code{ks::kcde}
#'       used in (Park and Gupta, 2012).}
#'   }
