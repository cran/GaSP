#' Predict from a \code{GaSPModel} object.
#' @inheritParams CrossValidate
#' @param x_pred A data frame containing the values of the input variables
#' at which to predict the output.
#' @param generate_coefficients A boolean indicating whether
#' coefficients for further external predictions are generated.
#' @return A list with the following elements:
#' \item{y_pred}{A data frame with two columns: the predictions
#' \code{Pred} and their standard errors \code{SE}.}
#' \item{pred_coeffs}{A vector of coefficients for further predictions;
#' \code{NULL} if \code{generate_coefficients} is \code{FALSE}.}
#'
#' @note The vector of prediction coefficients in \code{pred_coeffs}
#' can be used as follows.  Let \eqn{c} denote the coefficients and let
#' \eqn{r} denote a vector with element \eqn{i} containing the correlation
#' between the output at a given new point and the output at training point \eqn{i}.
#' Then the prediction for the output at the new point is the dot product
#' of \eqn{c} and \eqn{r}.
#'
#' \code{\link{RMSE}} computes the root mean squared error
#' of the predictions.
#' \code{\link{PlotPredictions}} and \code{\link{PlotResiduals}}
#' plot the predictions or their residuals;
#' \code{\link{PlotStdResiduals}} and \code{\link{PlotQQ}}
#' plot the standardized residuals.
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
#' }
#'
#' borehole_pred <- Predict(
#'   GaSP_model = borehole_fit,
#'   x_pred = borehole$x_pred,
#'   generate_coefficients = TRUE
#' )
#' @export
Predict <- function(GaSP_model,
                    x_pred,
                    generate_coefficients = c(FALSE, TRUE)) {
  if (class(GaSP_model) != "GaSPModel") {
    stop("Not a GaSP Model.", call. = FALSE)
  }
  if (!is.logical(generate_coefficients)) {
    stop("'generate_coefficients' must be logical.", call. = FALSE)
  }
  if (all(generate_coefficients == c(FALSE, TRUE))) {
    stop("'generate_coefficients' must be specified.", call. = FALSE)
  }
  backup_options <- options()
  on.exit(options(backup_options))
  Check <- .ErrorSummarizer()
  validate <- .GaSPModelValidate(GaSP_model, Check)
  x <- data.frame(lapply(validate$x, as.double))
  x_pred <- .XPredCheck(x_pred, x, Check)
  .ErrorOut(Check)
  y <- as.double(validate$y)
  cor_family <- .ArgtoInt(GaSP_model$cor_family)
  random_error <- validate$random_error
  reg_model <- GaSP_model$reg_model
  sp_model <- GaSP_model$sp_model
  sp_var <- GaSP_model$sp_var
  error_var <- GaSP_model$error_var
  cor_par <- validate$cor_par
  result <- .Call("predict", reg_model, sp_model, random_error, cor_family, x, y, x_pred, generate_coefficients, sp_var, error_var, cor_par)
  pred <- NULL
  pred$y_pred <- result[[1]]
  pred$pred_coeffs <- result[[2]]
  pred
}
