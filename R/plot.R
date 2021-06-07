#' Plot true versus predicted output.
#'
#' Plot true versus predicted output (response)
#' made by \code{Predict} or \code{CrossValidate}.
#'
#' @param y_pred A data frame of predicted output values
#' made by \code{Predict} or \code{CrossValidate}.
#' @param y A vector of length equal to the number of rows in \code{y_pred}
#' containing the true output values.
#' @param y_name An optional character string containing the output variable name
#' (for labels).
#' @param y_units An optional character string constaining the units
#' of the output variable (for labels).
#' @param title A character string for the name of the function generating
#' the predictions (for an appropriate title):
#'   "Predict" from \code{\link{Predict}} or
#'   "CrossValidate" from \code{\link{CrossValidate}};
#'   "" for no title.
#' @param pch Plotting symbol for \code{plot}; default is open circle.
#' @return No return value, generates plots.
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
#' PlotPredictions(borehole_cv, y, title = "CrossValidate")
#'
#' PlotPredictions(borehole_pred$y_pred, borehole$y_true, title = "Predict")
#' @export
PlotPredictions <- function(y_pred, y, y_name = "y", y_units = "",
                            title = c("Predict", "CrossValidate"), pch = 1) {
  backup_options <- options()
  on.exit(options(backup_options))
  results <- .PlotCheck(y, y_pred, y_name, y_units)
  y <- results$y
  y_pred <- results$y_pred
  y_name <- results$y_name
  y_units <- results$y_units
  if (y_units != "") {
    y_units <- paste("(", y_units, ")", sep = "")
  }
  x_label <- paste("Predicted", y_name, y_units, sep = " ")
  y_label <- paste(y_name, y_units)

  # Same scales for both axes.
  pred_col <- y_pred[[1]]
  neither.na <- !is.na(y) & !is.na(pred_col)

  # Remove NA's.
  pred_col <- pred_col[neither.na]
  y <- y[neither.na]

  xy.range <- range(y, pred_col)
  plot(pred_col, y,
    type = "p",
    xlab = x_label, ylab = y_label,
    xlim = xy.range, ylim = xy.range, pch = pch
  )

  # Add y = x line.
  abline(0, 1)

  # Add title
  if (title != "") {
    title <- .ArgtoInt(match.arg(title))
    if (title) {
      title(main = paste("True", y_name, "versus CV prediction"))
    } else {
      title(main = paste("True", y_name, "versus prediction"))
    }
  }
}

#' Plot residuals versus each input variable.
#' @inheritParams PlotPredictions
#' @param x A data frame with number of rows equal to the number of rows in
#' \code{y_pred} containing the input (explanatory) variables.
#' @param x_units An optional vector of character strings containing the units
#' of the input variables in \code{x} (for labels).
#' @return No return value, generates plots.
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
#' borehole_cv <- CrossValidate(borehole_fit)
#' }
#' PlotResiduals(x, borehole_cv, y)
#' @export
PlotResiduals <- function(x, y_pred, y, x_units = NULL, y_name = "y",
                          y_units = "", pch = 1) {
  backup_options <- options()
  on.exit(options(backup_options))
  results <- .PlotCheck(y, y_pred, y_name, y_units, x, x_units)
  y <- results$y
  y_pred <- results$y_pred
  y_name <- results$y_name
  y_units <- results$y_units
  if (y_units != "") {
    y_units <- paste("(", y_units, ")", sep = "")
  }
  x <- results$x
  x_units <- results$x_units
  if (!is.null(x_units)) {
    x_units <- paste("(", x_units, ")", sep = "")
  }
  x_name <- colnames(x)
  pred_col <- y_pred[[1]]
  res <- y - pred_col
  y_label <- paste(y_name, "residual", y_units, sep = " ")

  res_range <- range(res, na.rm = TRUE)
  res_na <- is.na(res)
  res[res_na] <- res_range[1] - 0.1 *
    (res_range[2] - res_range[1])
  res_range <- range(res)

  for (k in 1:ncol(x))
  {
    x_range <- range(x[[k]], na.rm = TRUE)
    plot(x[!res_na, k], res[!res_na],
      type = "p",
      xlab = paste(x_name[k], x_units[k]),
      ylab = y_label, xlim = x_range,
      ylim = res_range, pch = pch
    )
    if (length(res[res_na]) > 0) {
      # NAs to plot.
      points(x[res_na, k], res[res_na], pch = ".")
    }
    title(main = paste(y_name, "residual versus", x_name[k]))
  }
}

#' Normal Q-Q plot.
#'
#' Normal Q-Q plot of the standardized residuals
#' of predictions from \code{Predict} or \code{CrossValidate}.
#' @inheritParams PlotPredictions
#' @return No return value, generates plots.
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
#' borehole_cv <- CrossValidate(borehole_fit)
#' }
#' PlotQQ(borehole_cv, y)
#' @export
PlotQQ <- function(y_pred, y, y_name = "y") {
  backup_options <- options()
  on.exit(options(backup_options))
  results <- .PlotCheck(y, y_pred, y_name)
  y <- results$y
  y_pred <- results$y_pred
  y_name <- results$y_name
  pred_col <- y_pred[[1]]
  se_col <- y_pred[[2]]

  neither.na <- !is.na(y) & !is.na(pred_col) & !is.na(se_col)

  se_col <- se_col[neither.na]
  pred_col <- pred_col[neither.na]
  y <- y[neither.na]

  stand_res <- (y - pred_col) / se_col

  res_range <- range(stand_res, -3, 3)
  qqnorm(stand_res,
    xlab = "Standard normal quantile",
    ylab = paste(y_name, "standardized residual"),
    xlim = res_range, ylim = res_range
  )
  abline(0, 1)
}


#' Plot standardized residuals versus predictions.
#'
#' Plot standardized residuals versus predictions made by \code{Predict} or \code{CrossValidate}.
#' @inheritParams PlotPredictions
#' @return No return value, generates plots.
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
#' borehole_cv <- CrossValidate(borehole_fit)
#' }
#' PlotStdResiduals(borehole_cv, y, title = "CrossValidate")
#' @export
PlotStdResiduals <- function(y_pred, y, y_name = "y", y_units = "",
                             title = c("Predict", "CrossValidate"), pch = 1) {
  backup_options <- options()
  on.exit(options(backup_options))
  results <- .PlotCheck(y, y_pred, y_name, y_units)
  y <- results$y
  y_pred <- results$y_pred
  y_name <- results$y_name
  y_units <- results$y_units
  if (y_units != "") {
    y_units <- paste("(", y_units, ")", sep = "")
  }
  x_label <- paste("Predicted", y_name, y_units, sep = " ")

  pred_col <- y_pred[[1]]
  se_col <- y_pred[[2]]

  neither.na <- !is.na(y) & !is.na(pred_col) & !is.na(se_col)

  # Remove NA's.
  se_col <- se_col[neither.na]
  pred_col <- pred_col[neither.na]
  y <- y[neither.na]

  stand_res <- (y - pred_col) / se_col

  plot(pred_col, stand_res,
    type = "p", ylim = range(-4, 4, stand_res),
    xlab = x_label, ylab = "Standardized residual", pch = pch
  )

  # Add lines at +/-2 and +/-3
  abline(-2, 0, lty = 3, col = "red")
  abline(2, 0, lty = 3, col = "red")
  abline(-3, 0, lty = 2, col = "red")
  abline(3, 0, lty = 2, col = "red")

  if (title != "") {
    title <- .ArgtoInt(match.arg(title))
    if (title) {
      title(main = paste(y_name, "standardized residual versus CV prediction"))
    } else {
      title(main = paste(y_name, "standardized residual versus prediction"))
    }
  }
}
