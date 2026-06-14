#' @importFrom stats var ks.test
#' @importFrom nortest ad.test
#' @importFrom extRC dfm
#'
# Forward Orthogonal Deviations (FOD) and GLS transformation

# Applying combined FOD and GLS to a single panel's data vector

copulaPanelFod <- function(x){ #GLS corrected first-differences with a spherical error
  #covariance
  n <- length(x)
  if (n == 1L) return (NA_real_)
  D <- extRC::dfm(n) #first-diffrence matrix dim (n-1) x n
  A <- chol(solve(D %*% t(D))) %*% D #GLS
  c(as.vector(A %*% as.numeric(x)), NA_real_) #applying to x and append NA
}

#Applying the FOD & GLS transformation panel by panel to all columns

copulaPanelTransform <- function(data, col.panel, col.time, col.var){

  panel <- split(data, data[[col.panel]])

  resultlist <- lapply(panel, function(panel_i){
    panel_i <- panel_i[order(panel_i[[col.time]]), , drop = FALSE] #Within each panel, we order by time
    #this is needed for the FD to be well-defined
    n_i <- nrow(panel_i)
    if (n_i <= 1L){
      return(NULL)
    }

    #applying copulaPanelFod column by column

    matrix <- as.matrix(panel_i[, col.vars, drop = FALSE])
    out <- as.data.frame(apply(matrix, 2, copulaPanelFod))
    colnames(out) <- col.vars

    #keeping the panel IDs for bootstrap
    out[[col.panel]] <- panel_i[[col.panel]][1L]
    out[stats::complete.cases(out), , drop = FALSE]#removing NA row
  })

  do.call(rbind, Filter(Negate(is.null), resultlist))

}

# Identification checks

#1. Jarque-bera test
#using the Jarque-Bera test on the structural residuals.
#if the P-value is low, there may be error distribution issues. It means non-symmetrical
# residuals

copulaPanel_JB <- function(residuals){
  tsoutliers::JarqueBera.test(residuals)[[2L]]$p.value
}

#2. the Anderson-Darling test for each endo regressor
# if there is high p-value, it is close to normal. This can lead to identification issues.

copulaPanel_AD <- function(P){
  apply(P, 2, function(x) suppressWarnings(nortest::ad.test(scale(x))$p.value))
}

#3. KS test between the scaled residuals and each endo regressor
#if the p-value is high, it means the distribution are the same. This can lead to identification
#issues

copulaPanel_KS <- function(residuals, P){
  apply(P, 2, function(x){
    suppressWarnings(stats::ks.test(scale(residuals), scale(x))$p.value)
  })
}

#BOOTSTRAP
#pair bootstrap approach (from eq. 26)
#cross sectional pairs bootstrap

bootstrapPanel <- function(fn.fit, data.fod, col.panel, num.boots, coef.names, verbose){

  panel.id <- unique(data.fod[[col.panel]])
  panel.length <- length(panel.id)
  boots.params <- matrix (NA_real_, nrow = num.boots, ncol= length(coef.names))
  colnames(boots.params) <- coef.names
  n.failed <- 0L

  if (verbose) message("Running", num.boots, "bootstrap replications...")

  for (b in seq_len(num.boots)){
    sample.id <- sample(panel.id, panel.length, replace =TRUE)
    data.boot <- do.call(rbind, lapply(sample.id, function(id){
      data.fod[data.fod[[col.panel]] == id, ,drop = FALSE]
    }))

    result <- tryCatch(
      supressWarnings(fn.fit(data.boot)),
      error = function(e) NULL
    )
    if (is.null(result) || anyNA(result)){
      n.failed <- n.failed + 1L
    } else{
      boots.params[b, ] <- result
    }
  }

  if (verbose && n.failed > 0L){
    message(n.failed, "of", num.boots, "bootstrap replications failed (", round(100*n.failed/num.boots, 1), "%)
    and were discarded")
  }

  list(boots.params = boots.params,
  n.attempted = num.boots,
  n.failed = n.failed
  )
}


