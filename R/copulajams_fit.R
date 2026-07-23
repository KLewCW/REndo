#' @importFrom Formula as.Formula
#' @importFrom stats lm model.frame model.matrix terms formula update reformulate
copulajams_fit <- function(F.formula, data, labels.endo, labels.exo, cdf) {
  # Build a single model.frame/matrix from which all parts (W,P,Z) are read from to
  # guarantee they share row and col order
  mf.main <- model.frame(F.formula, lhs = 1, rhs = 1, data = data, na.action = na.fail)
  # get terms from mf to also contain dataClasses
  tr.main <- terms(mf.main)
  X.main <- model.matrix(tr.main, data = mf.main)

  # get term labels in correct order for the "assign" mapping of the model.matrix
  # (cannot use labels.main or union(labels.endo, labels.exo))
  term.labels <- labels(tr.main)
  # for each column in model matrix: Index of the term.label that generated it
  mm.assign <- attr(X.main, "assign")

  # find exo labels which are factors or not
  is.exo.factor <- vapply(
    labels.exo,
    FUN = function(lbl.exo) {
      dc <- attr(tr.main, "dataClasses")
      return(dc[label_to_colname(lbl.exo)] %in% c("factor", "ordered"))
    },
    FUN.VALUE = logical(1)
  )

  labels.exo.factor <- labels.exo[is.exo.factor]
  labels.exo.cont <- labels.exo[!is.exo.factor]

  X.endo <- X.main[,
    mm.assign %in% match(labels.endo, table = term.labels),
    drop = FALSE
  ]
  X.exo.cont <- X.main[,
    mm.assign %in% match(labels.exo.cont, table = term.labels),
    drop = FALSE
  ]

  # message("X.main")
  # str(X.main)
  # message("term.labels", toString(term.labels))
  # message("labels.endo", toString(labels.endo))
  # message("mm.assign")
  # str(mm.assign)
  # message("X.endo")
  # str(X.endo)
  # message("labels.exo.cont: ", labels.exo.cont)
  # message("X.exo.cont")
  # str(X.exo.cont)

  # if (nrow(X.endo) < length(labels.endo)) {
  #   stop(
  #     "Bootstrap sample dropped at least one endogenous regressor.",
  #     "This can happen when a regressor becomes constant in a resample.",
  #     call. = FALSE
  #   )
  # }

  # Checking if there is any factor variables among exogenous regressors
  # if the factor Z is present, then it requires stratified correction per level (see equation 20 and 21)
  # if there is no factor Z, variance-covariance correction (eq. 17 to 19)
  if (length(labels.exo.factor) == 0) {
    # case no Z
    cop.terms <- copulajams_correction_continuous(
      P = X.endo,
      W = X.exo.cont,
      cdf = cdf
    )
  } else {
    # case where Z is present
    # Z is only used to apply-by-factor: Has to remain factors and not made to dummies
    colnames.exo.factors <- vapply(
      labels.exo.factor,
      FUN = label_to_colname,
      FUN.VALUE = character(1)
    )
    df.Z <- mf.main[, colnames.exo.factors, drop = FALSE]
    colnames(df.Z) <- labels.exo.factor

    cop.terms <- copulajams_correction_discrete(
      P = X.endo,
      W = X.exo.cont,
      df.Z = df.Z,
      cdf = cdf
    )
  }

  # 2nd stage: OLS ---------------------------------------------------------------------

  res.2nd.stage <- copula_fit_2ndstage(
    F.formula = F.formula,
    data = data,
    cop.terms = cop.terms
  )

  return(list(
    res.augmented = res.2nd.stage$res.augmented,
    labels.pcop = res.2nd.stage$labels.pcop
  ))
}


#' @importFrom stats qnorm cov
copulajams_correction_continuous <- function(P, W, cdf) {
  # cols must be named and may not contain intercept
  stopifnot(!is.null(colnames(P)) & !("" %in% colnames(P)))
  stopifnot(!is.null(colnames(W)) & !("" %in% colnames(W)))
  stopifnot(!"(Intercept)" %in% colnames(P))
  stopifnot(!"(Intercept)" %in% colnames(W))

  #computing the copula correction terms - continuous exogenous case
  # equation 17 : C(P_i, W_i), sigma_(C(P), C(W)} is the variance-covariance matrix
  #of (C(P), W(P)) (not the correlation matrix)
  # the final (I_{dp}, 0)' projection keeps only the d_p columns which correspond to the endo
  #regressors and giving one copula term per endo regressors

  PW <- cbind(P, W)

  #step 1: applying CDF to all regressors (P&W)
  P.star <- copula_pstar(P = PW, cdf = cdf)

  #step 2: applying qnorm to obtain normal scores C(P) and C(W)
  C.PW <- apply(P.star, 2, qnorm) #using equation C(P_i) = Phi^{-1}(F_hat(P_i)) for each regressor

  #step 3: Estimate variance-covariance matrix \hat{sigma}_{C(P), C(W)}
  Sigma.hat <- cov(C.PW) #equation 21: mentioned that it is the variance-covariance so cov should be used
  #instead of cor(). i think ?

  #step 4: inverse of the variance-covariance matrix:
  Sigma.inverse <- solve(Sigma.hat)

  #step 4:finding the d_P copula term to project onto endogenous regressor columns
  # (I_{dp}, 0_{dw x dp})' selects the d_P endo columns of Sigma^{-1}
  # giving the d_P copula terms (one per endo reg)
  endo.index <- seq_len(ncol(P)) # P is in the first ncol() columns

  #C(P_i, W_i) = (C(P_i)', C(W_i)') \hat(Sigma^{-1}_{C(P), C(W)}) (I_{dp}, 0_{dw x dp})'
  P.cop <- C.PW %*% Sigma.inverse[, endo.index, drop = FALSE]
  P.cop <- as.matrix(P.cop)
  colnames(P.cop) <- paste0(colnames(P), "_cop")

  return(P.cop)
}


#' @importFrom stats cov qnorm
copulajams_correction_discrete <- function(
  P,
  W,
  df.Z, # has to be factors, not dummies
  cdf
) {
  stopifnot(!is.null(colnames(df.Z)))
  l.results <- list()

  # For each factor separately: For each level: Process data
  for (label.Z.i in colnames(df.Z)) {
    # Data of single factor (drop: need vector not df to read with levels())
    Z.i <- df.Z[, label.Z.i, drop = TRUE]

    # levels() also applies for ordered factors
    for (lvl.i in levels(Z.i)) {
      # find where data with this level is
      idx.rows <- Z.i == lvl.i

      # check num obs before doing expensive data subset
      if (sum(idx.rows) <= 3) {
        # message("Has less than 3 obs")
        # 3 is number checked by Haschka's repo line 190.
        #did not find paper backing this up ?
        next
      }

      P.sub <- P[idx.rows, , drop = FALSE]
      W.sub <- W[idx.rows, , drop = FALSE]

      # message("P")
      # str(P)
      # message("W")
      # str(W)
      # message("P.sub")
      # str(P.sub)
      # message("W.sub")
      # str(W.sub)
      # message("has variation")
      # str(apply(P.sub, 2, function(col) length(unique(col)) > 1))

      #skip if too few observations or no variation in endogenous regressors
      # has.variation <- any(sapply(
      #   subdat1[, names.endo.regs, drop = FALSE],
      #   function(col) length(unique(col)) > 1
      # ))
      # if (nrow(subdat1) <= 3 || !has.variation) {
      has.variation <- any(apply(P.sub, 2, function(col) length(unique(col)) > 1))
      if (!has.variation) {
        # message("Has no variation")
        next
      }

      ##usually need at least p+1 observations to estimate a pxp cov matrix
      # more generally, instead of '3', we could have tried the most conservative minimum of
      # max(10, 2*p) for reliable estimation
      #min.observation.needed <- max(10L, 2L * length(cols.use))
      #if(nrow(subdat1) <= min.observation.needed || !has.variation){
      # then maybe issue a warning here... saying that the factor level 'lvl' of the variable 'var'
      # does not have enough observation for a reliable covariance estimation and we are skipping this level.

      #now use non-factor columns only for the CDF transformation.
      #factor vairables cannot enter the CDF transformation as they are discrete with no meaningful continuous CDF

      # cols.use <- setdiff(c(names.endo.regs, names.exo.regs), factor.vars)
      # subdat2 <- as.matrix(subdat1[, cols.use, drop = FALSE])
      #
      # subdat2 <- subdat2[,
      #                    apply(subdat2, 2, function(x) length(unique(x)) > 1),
      #                    drop = FALSE
      # ] #keeping only columns with variation with subset
      #
      # if (ncol(subdat2) == 0) {
      #   next
      # }

      #Using the steps from equation 21 again:

      # Other than for the continuous-only case, we do not want an error (like
      # non-invertible matrix) in a single strata to derail the whole pcop generation
      # process.
      #
      # So far, continue with all 0 only for errors related to matrix inversion.
      # Let all other errors & warnings go through
      # TODO: maybe wrap solve() itself + throw custom error to identify more robustly
      cop.terms.sub <- tryCatch(
        copulajams_correction_continuous(P = P.sub, W = W.sub, cdf = cdf),
        error = function(e) {
          # non-invertible signals an error containing "singular"
          if(grepl("singular", conditionMessage(e), ignore.case = TRUE)){
            return(NULL)
          }
          # fail / "re-throw" any other error
          stop(e)
        }
      )

      # there was an issue (with inverting the matrix)
      if (is.null(cop.terms.sub)) {
        #return zero correction terms
        # create placeholder
        cop.terms.sub <- matrix(0, nrow = nrow(P.sub), ncol = ncol(P.sub))

      }

      # naming correction terms with factor level info
      # K.actual <- ncol(P.cop)
      # endo.for.names <- names.endo.regs[seq_len(K.actual)]
      # colnames(P.cop) <- paste0(endo.for.names, "_", var, "_", lvl, "_cop")

      # Expanding back to a full dataset
      # I(Z_i = z) from eq. 20. zero attributed for observations not in this level


      # As many rows as full data (P) but only as many cols as correction applied
      cop.terms.full <- matrix(0, nrow = nrow(P), ncol = ncol(cop.terms.sub))
      cop.terms.full[idx.rows, ] <- cop.terms.sub


      # name columns according to strata (name)
      name.strata <- paste(label.Z.i, lvl.i, sep = "_")
      colnames(cop.terms.full) <- make.names(paste0(name.strata, "_", colnames(cop.terms.sub)))
      l.results[[name.strata]] <- cop.terms.full
    }
  }

  # str(l.results)

  if (length(l.results) == 0) {
    stop(
      "No valid factor-level subsets found for correction term estimation. ",
      "Check that factor exogenous regressors have sufficient observations ",
      "per level (> 3) and variation in the endogenous regressors.",
      call. = FALSE
    )
  }

  cop.terms.results <- do.call(cbind, l.results)
  # message("cop.terms.results")
  # str(cop.terms.results)
  return(cop.terms.results)
}
