#' @importFrom utils getS3method
checkinput_copula2scopenp_npcdistbwargs <- function(npcdistbw.args) {
  # No checks if parameter default
  if (identical(npcdistbw.args, list())) {
    return(c())
  }

  err.msg <- c()

  # Plain list
  if (!identical(class(npcdistbw.args), "list")) {
    return("Parameter `npcdistbw.args` must be a plain list.")
  }

  # Must be named
  nms <- names(npcdistbw.args)
  if (is.null(nms) || any(nms == "")) {
    return("All elements of 'npcdistbw.args' must be named.")
  }

  # Reserved
  reserved.nms <- intersect(nms, c("xdat", "ydat"))
  if (length(reserved.nms) > 0L) {
    err.msg <- c(
      err.msg,
      paste0(
        "`npcdistbw.args` may not contain reserved argument(s): ",
        paste(reserved.nms, collapse = ", ")
      )
    )
  }

  valid.nms <- names(formals(getS3method(f = "npcdistbw", class = "default")))
  unknown.nms <- setdiff(nms, valid.nms)
  if (length(unknown.nms) > 0L) {
    unknown.str <- paste(unknown.nms, collapse = ", ")
    warning(
      "Unknown argument(s) in `npcdistbw.args`: ",
      unknown.str,
      call. = FALSE,
      immediate. = TRUE
    )
  }

  return(err.msg)
}
