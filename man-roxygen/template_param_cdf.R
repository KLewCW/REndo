#' @param cdf Character string specifying the method used to estimate the
#'   marginal distribution functions of the endogenous regressors.
#'   One of the following:
#' \describe{
#'     \item{\code{"adj.ecdf"}}{Adjusted empirical CDF with midrank correction
#'       (Liengaard et al., 2024).
#'       \deqn{\hat{F}(x_i) = \frac{1}{2n} + \frac{n-1}{n^2} \sum_{i=1}^{n} \mathbf{1}(P_i \leq x)}
#'       where n is the sample size, \eqn{P_i} is the endogenous regressor for observation i = 1,..,n and \eqn{\mathbf{1} (.)} is the
#'       indicator function.
#'       This is a sample size-dependent adjustment of the ECDF. Correction term is acquired by minimising
#'       the true CDF's MSE, using the standard ECDF as a predictor in a linear regression.
#'       Keeps values strictly inside \eqn{(0,1)} but retains the ECDF's desirable large sample properties.}
#'     \item{\code{"resc.ecdf"}}{Rescaled empirical CDF via
#'       \code{copula::pobs}. \deqn{\hat{F}(x_i) = \frac{\text{rank}(x_i)}{n + 1}}
#'        Keeps values strictly inside \eqn{(0,1)} without an arbitrary boundary constant.
#'       Breitung et al. (2024) adopts this as their theoretical recommendation.}
#'      \item{\code{"ecdf"}}{<%= if (exists("ecdf.text")) ecdf.text else "Empirical CDF with boundary replacement (Becker et al., 2022)." %>}
#'     \item{\code{"kde"}}{Integral of a density estimator via \code{ks::kcde}
#'       used in (Park and Gupta, 2012).}
#'   }
