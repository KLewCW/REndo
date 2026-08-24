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
    return(paste0(
      "Parameters 'num.iterations' (", num.iterations, "), 'burnin' (", burnin,
      ") and 'thin' (", thin, ") retain only ", num.remaining,
      " draw(s). At least 2 are required."))
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
