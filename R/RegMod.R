# Function to read the usual R model description and stores them in a data frame
.RegMod <- function(reg_model) {
  rterms <- terms(reg_model)
  intercept <- toString(attr(rterms, "intercept"))
  if (intercept == "0") {
    intercept <- NULL
  }
  rmod <- attr(rterms, "term.labels")

  if (length(rmod) == 0) {
    rmod <- NULL
  } else {
    rmod <- gsub("I(", "", rmod, fixed = T)
    rmod <- gsub(")", "", rmod, fixed = T)
  }

  terms <- c(intercept, rmod)

  m <- data.frame(Terms = terms, stringsAsFactors = FALSE)
  m
}

.SPModHelper <- function(sp_model, x) {
  sp <- NULL
  if (is.null(sp_model)) {
    sp_model_df <- NULL
    x_colname <- colnames(x)
    sp_model_df <- data.frame(Terms = x_colname, stringsAsFactors = FALSE)
    sp <- sp_model_df
  } else {
    sp <- .SPMod(sp_model)
  }
  return(sp)
}

.SPMod <- function(sp_model) {
  rterms <- terms(sp_model)
  rmod <- attr(rterms, "term.labels")

  if (length(rmod) == 0) {
    rmod <- NULL
  }
  else {
    rmod <- gsub("I(", "", rmod, fixed = T)
    rmod <- gsub(")", "", rmod, fixed = T)
  }

  terms <- c(rmod)

  m <- data.frame(Terms = terms, stringsAsFactors = FALSE)
  m
}
