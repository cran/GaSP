#' Calculate the root mean squared error (RMSE) of prediction
#' @param y_pred A vector of predicted output values.
#' @param y_true A vector of true output values.
#' @param normalized An optional boolean: if \code{TRUE}, the RMSE is
#' normalized by dividing it by the standard deviation of \code{y_true}.
#' @return The RMSE or normalized RMSE.
#' @examples
#' \dontshow{
#' x <- borehole$x
#' y <- borehole$y
#' theta <- c(
#'   5.767699e+01, 0.000000e+00, 0.000000e+00, 1.433571e-06,
#'   0.000000e+00, 2.366557e-06, 1.695619e-07, 2.454376e-09
#' )
#' alpha <- c(
#'   1.110223e-16, 0.000000e+00, 0.000000e+00, 0.000000e+00,
#'   0.000000e+00, 0.000000e+00, 2.494862e-03, 0.000000e+00
#' )
#' cor_par <- data.frame(Theta = theta, Alpha = alpha)
#' rownames(cor_par) <- colnames(borehole$x)
#' sp_var <- 38783.7
#' borehole_fit <- GaSPModel(
#'   x = borehole$x, y = borehole$y,
#'   reg_model = ~1, cor_family = "PowerExponential",
#'   cor_par = cor_par, random_error = FALSE,
#'   sp_var = sp_var
#' )
#' borehole_pred <- Predict(
#'   GaSP_model = borehole_fit,
#'   x_pred = borehole$x_pred,
#'   generate_coefficients = TRUE
#' )
#' borehole_cv <- CrossValidate(borehole_fit)
#' }
#' RMSE(borehole_pred$y_pred$Pred, borehole$y_true)
#'
#' RMSE(borehole_cv$Pred, y)
#' @export
RMSE <- function(y_pred, y_true, normalized = FALSE) {
  backup_options <- options()
  on.exit(options(backup_options))
  v <- .YPredCheck(y_pred, y_true, normalized)
  y_pred <- v$y_pred
  y_true <- v$y_true
  rmse <- sqrt(sum((y_true - y_pred)^2, na.rm = TRUE) / sum(!is.na(y_true)))
  if (normalized) {
    rmse.mean <- sqrt(sum((y_true - mean(y_true))^2, na.rm = TRUE) / sum(!is.na(y_true)))
    rmse <- rmse / rmse.mean
  }
  return(rmse)
}
