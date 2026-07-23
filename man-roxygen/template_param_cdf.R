#' @param cdf Character string specifying the method used to estimate the
#'   marginal distribution functions of the endogenous regressors.
#'   See details section.
#' \describe{
#'    \item{\code{"adj.ecdf"}}{Adjusted empirical CDF with midrank correction (Liengaard et al., 2024)}
#'    \item{\code{"resc.ecdf"}}{Rescaled empirical CDF=rank(x)/(n+1). <%= if (exists("resc.ecdf.addition")) resc.ecdf.addition else "" %>}
#'    \item{\code{"ecdf"}}{Empirical CDF with boundary replacement (Becker et al., 2022).}
#'    \item{\code{"kde"}}{Integral of a density estimator via \code{ks::kcde} used in (Park and Gupta, 2012).}
#'   }
#'
