.ErrorSummarizer <- function() {
  backup_options <- options()
  on.exit(options(backup_options))
  check <- new.env()
  assign("error", NULL, envir = check)
  options(warn = 1)
  options(error = NULL)
  check
}

.ErrorOut <- function(check) {
  error_l <- length(check$error)
  if (error_l > 0) {
    stop(paste0(c("", paste0(
      1:error_l,
      ": ", check$error
    )), sep = "\n"), call. = FALSE)
  }
}

.Error <- function(err_msg, check) {
  check$error <- c(check$error, err_msg)
}
