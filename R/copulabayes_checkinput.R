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


checkinput_copulabayes_method <- function(method, allowed.methods) {
  # Default value of the param is the full vector of allowed values (match.arg() then
  # selects the first)
  if (identical(method, allowed.methods)) {
    return(c())
  }

  if (!is.character(method) || length(method) != 1) {
    return("Parameter 'method' must be a single character string.")
  }

  if (!(method %in% allowed.methods)) {
    return(paste0("Parameter 'method' must be one of: ", toString(allowed.methods), "."))
  }

  return(c())
}
