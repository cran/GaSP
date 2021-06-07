.ArgtoInt <- function(x) {
  if (is.character(x) & length(x) == 1) {
    if (all(x == "Matern") || all(x == "Posterior") ||
      all(x == "CV") || all(x == "CrossValidate")) {
      return(1)
    }
    if (all(x == "PowerExponential") || all(x == "Likelihood") ||
      all(x == "Objective") || all(x == "Predict")) {
      return(0)
    }
  }
}
