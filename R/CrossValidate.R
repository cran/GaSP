#' Cross-validated predictions for a \code{GaSPModel} object.
#'
#' Compute leave-one-out cross-validated predictions for a \code{GaSPModel} object.
#'
#' @param GaSP_model Object of class \code{\link{GaSPModel}}.
#' @return A data frame with two columns: the cross-validated predictions
#' \code{Pred} and their standard errors \code{SE}.
#' @note \code{\link{RMSE}} computes the root mean squared error
#' of the predictions.
#' \code{\link{PlotPredictions}} and \code{\link{PlotResiduals}}
#' plot the predictions or their residuals;
#' \code{\link{PlotStdResiduals}} and \code{\link{PlotQQ}}
#' plot the stanadardized residuals.
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
#' borehole_cv <- CrossValidate(borehole_fit)
#' @export
CrossValidate <- function(GaSP_model) {
  if (!inherits(GaSP_model, 'GaSPModel')) {
    stop("Not a GaSP Model.")
  }
  backup_options <- options()
  on.exit(options(backup_options))
  Check <- .ErrorSummarizer()
  validate <- .GaSPModelValidate(GaSP_model, Check)
  x <- data.frame(lapply(validate$x, as.double))
  .ErrorOut(Check)
  y <- as.double(validate$y)
  cor_family <- .ArgtoInt(GaSP_model$cor_family)
  random_error <- validate$random_error
  reg_model <- GaSP_model$reg_model
  sp_model <- GaSP_model$sp_model
  sp_var <- GaSP_model$sp_var
  error_var <- GaSP_model$error_var
  cor_par <- validate$cor_par
  result <- .Call("crossvalidate", reg_model, sp_model, random_error, cor_family, x, y, sp_var, error_var, cor_par)
  result
}
