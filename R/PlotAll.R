#' Execute \code{PlotPredictions}, \code{PlotResiduals}, \code{PlotStdResiduals}, \code{PlotMainEffects}, and \code{PlotJointEffects}.
#'
#' Execute \code{PlotPredictions}, \code{PlotResiduals} and \code{PlotStdResiduals} (all applied to CV only),
#' \code{PlotMainEffects}, and \code{PlotJointEffects}.
#'
#' @param GaSP_model Object of class \code{\link{GaSPModel}},
#' the entire model will be verified but only \code{x} and \code{y} will be used.
#' @param cross_validation A data frame returned by \code{\link{CrossValidate}}.
#' @param visualization A list object returned by \code{\link{Visualize}}.
#' @return No return value, generates plots.
#' @inheritParams PlotJointEffects
#' @inheritParams PlotPredictions
#' @importFrom stats qnorm qqnorm qt spline terms
#' @importFrom graphics abline contour lines par points title plot
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
#' borehole_x_names <- colnames(x)
#' borehole_min <- c(0.05, 100.00, 63070.00, 990.00, 63.10, 700.00, 1120.00, 9855.00)
#' borehole_max <- c(0.15, 50000.00, 115600.00, 1110.00, 116.00, 820.00, 1680.00, 12045.00)
#' borehole_x_desc <- DescribeX(borehole_x_names, borehole_min, borehole_max)
#' borehole_vis <- Visualize(borehole_fit, borehole_x_desc)
#' borehole_cv <- CrossValidate(borehole_fit)
#' }
#' PlotAll(borehole_fit, borehole_cv, borehole_vis)
#' @export
PlotAll <- function(GaSP_model, cross_validation,
                    visualization,
                    y_name = "y", y_units = "", x_units = NULL, se_plot = TRUE,
                    y_values = NULL, se_values = NULL, pch = 1) {
  if (class(GaSP_model) != "GaSPModel") {
    stop("Not a GaSP Model.")
  }
  backup_options <- options()
  on.exit(options(backup_options))
  Check <- .ErrorSummarizer()
  .GaSPModelValidate(GaSP_model, Check)
  .ErrorOut(Check)

  x <- GaSP_model$x
  y <- GaSP_model$y
  anova_percent <- visualization$anova_percent
  main_effect <- visualization$main_effect
  joint_effect <- visualization$joint_effect

  PlotPredictions(cross_validation, y, y_name = y_name, y_units = y_units, title = "CrossValidate", pch = pch)

  PlotResiduals(x, cross_validation, y,
    y_name = y_name, y_units = y_units,
    x_units = x_units, pch = pch
  )

  PlotStdResiduals(cross_validation, y, y_name = y_name, y_units = y_units, title = "CrossValidate", pch = pch)

  PlotMainEffects(main_effect, anova_percent, y_name = y_name, y_units = y_units, x_units = x_units)

  PlotJointEffects(joint_effect, anova_percent,
    y_name = y_name, y_units = y_units,
    x_units = x_units, se_plot = se_plot,
    y_values = y_values, se_values = se_values
  )
}
