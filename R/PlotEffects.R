#' Plot the estimated main effects.
#'
#' @param main_effect A data frame from \code{\link{Visualize}}
#' with plotting coordinates for the estimated main effects.
#' @param anova_percent A data frame from \code{\link{Visualize} of ANOVA percentages}.
#' @param x_units An optional vector of character strings containing the units
#' of the input variables (for labels).
#' @param y_name An optional character string containing the output variable name
#' (for labels).
#' @param y_units An optional character string containing the units
#' of the output variable (for labels).
#' @return No return value, generates plots.
#' @details
#' Plots are sent to the active device.
#' Each plot shows an estimated main effect (red solid line)
#' and pointwise approximate 95% confidence limits (green dashed line).
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
#' borehole_x_names <- colnames(borehole$x)
#' borehole_min <- c(0.05, 100.00, 63070.00, 990.00, 63.10, 700.00, 1120.00, 9855.00)
#' borehole_max <- c(0.15, 50000.00, 115600.00, 1110.00, 116.00, 820.00, 1680.00, 12045.00)
#' borehole_x_desc <- DescribeX(borehole_x_names, borehole_min, borehole_max)
#' borehole_vis <- Visualize(borehole_fit, borehole_x_desc)
#' }
#' PlotMainEffects(borehole_vis$main_effect, borehole_vis$anova_percent)
#' @export
PlotMainEffects <- function(main_effect, anova_percent, x_units = NULL,
                            y_name = "y", y_units = "") {
  backup_options <- options()
  on.exit(options(backup_options))
  me_list <- .MainEffectExtract(main_effect)
  x <- me_list$x
  y <- me_list$y
  se <- me_list$se
  x_name <- me_list$x_name
  index_start <- me_list$index_start
  index_finish <- me_list$index_finish
  n_plots <- length(index_start)
  x_names <- .XNameExtract(anova_percent)
  check <- .ErrorSummarizer()
  check_result <- .EffectPlotCheck(
    anova_percent, main_effect, NULL,
    y_name, y_units, x_units, x_names, check
  )
  .ErrorOut(check)
  anova_percent <- check_result$anova_percent
  main_effect <- check_result$main_effect
  x_units <- check_result$x_units
  if (y_units != "") {
    y_units <- paste("(", y_units, ")", sep = "")
  }
  if (!is.null(x_units)) {
    x_units <- paste("(", x_units, ")", sep = "")
  }
  y_range <- range(
    y - 1.96 * se,
    y + 1.96 * se
  )
  for (i in 1:n_plots)
  {
    plot_cases <- index_start[i]:index_finish[i]

    plot(spline(x[plot_cases], y[plot_cases]),
      type = "l", lty = 1, col = "red",
      xlab = paste(x_name[i], x_units[.XUnitsFind(x_name[i], x_names)]),
      # xaxs = "e",
      ylab = paste(y_name, y_units),
      ylim = y_range
    )
    lines(spline(x[plot_cases], y[plot_cases] - 1.96 * se[plot_cases]),
      lty = 2, col = "green4"
    )
    lines(spline(x[plot_cases], y[plot_cases] + 1.96 * se[plot_cases]),
      lty = 2, col = "green4"
    )

    yx.title <- paste(y_name, "(", x_name[i], ") : ",
      .ANOVAFormat(anova_percent[x_name[i], "y"]),
      "%",
      sep = ""
    )
    title(yx.title)
  }
}

.MainEffectExtract <- function(main_effect) {
  x_name <- NULL
  y_name <- NULL
  index_start <- NULL
  index_finish <- NULL

  i <- 1
  while (i <= nrow(main_effect)) {
    x_name_this_plot <- as.character(main_effect[i, 1])
    j <- i
    while (j < nrow(main_effect) &&
      main_effect[j + 1, 1] == x_name_this_plot) {
      j <- j + 1
    }

    # The plot of y_name.this.plot versus x_name_this_plot
    # uses plotting coordinates in rows i through j.

    x_name <- c(x_name, x_name_this_plot)
    # y_name <- c(y_name, y_name.this.plot)
    index_start <- c(index_start, i)
    index_finish <- c(index_finish, j)

    i <- j + 1
  }

  return(list(
    x = as.numeric(main_effect[, 2]),
    y = as.numeric(main_effect[, 3]),
    se = as.numeric(main_effect[, 4]),
    x_name = x_name,
    index_start = index_start,
    index_finish = index_finish
  ))
}

#' Plot the estimated joint effects.
#' @inheritParams PlotMainEffects
#' @param joint_effect A data frame from \code{\link{Visualize}}
#' with plotting coordinates for the estimated joint effects.
#' @param se_plot An optional boolean indicating whether to make
#' standard-error contour plots.
#' @param y_values An optional vector of contour values for the estimated
#' joint effects.
#' @param se_values An optional vector of contour values for the standard
#' errors.
#' @return No return value, generates plots.
#' @details
#' Plots are sent to the active device.
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
#' borehole_x_names <- colnames(borehole$x)
#' borehole_min <- c(0.05, 100.00, 63070.00, 990.00, 63.10, 700.00, 1120.00, 9855.00)
#' borehole_max <- c(0.15, 50000.00, 115600.00, 1110.00, 116.00, 820.00, 1680.00, 12045.00)
#' borehole_x_desc <- DescribeX(borehole_x_names, borehole_min, borehole_max)
#' borehole_vis <- Visualize(borehole_fit, borehole_x_desc)
#' }
#' PlotJointEffects(borehole_vis$joint_effect, borehole_vis$anova_percent)
#' @export
PlotJointEffects <- function(joint_effect, anova_percent, x_units = NULL,
                             y_name = "y", y_units = "",
                             se_plot = TRUE,
                             y_values = NULL, se_values = NULL) {
  backup_options <- options()
  on.exit(options(backup_options))
  check <- .ErrorSummarizer()
  x_names <- .XNameExtract(anova_percent)
  check_result <- .EffectPlotCheck(anova_percent, NULL, joint_effect, y_name, y_units, x_units, x_names, check)
  check_result_add <- .JointEffectAdditCheck(se_plot, y_values, se_values)
  .ErrorOut(check)
  anova_percent <- check_result$anova_percent
  joint_effect <- check_result$joint_effect
  x_units <- check_result$x_units
  se_plot <- check_result_add$se_plot
  y_values <- check_result_add$y_values
  se_values <- check_result_add$se_values
  if (y_units != "") {
    y_units <- paste("(", y_units, ")", sep = "")
  }
  i <- 1
  while (i <= nrow(joint_effect)) {
    x1_name <- as.character(joint_effect[i, 1])
    x2_name <- as.character(joint_effect[i, 2])
    y <- as.numeric(joint_effect[, 5])
    se <- as.numeric(joint_effect[, 6])

    j <- i
    while (j < nrow(joint_effect) &&
      joint_effect[j + 1, 1] == x1_name &&
      joint_effect[j + 1, 2] == x2_name) {
      j <- j + 1
    }

    # The plot of y variable y_name versus x1_name and
    # x2_name uses plotting coordinates in rows i through j.
    x1 <- as.numeric(joint_effect[i:j, 3])
    x2 <- as.numeric(joint_effect[i:j, 4])
    x1 <- sort(unique(x1))
    x2 <- sort(unique(x2))

    y <- matrix(y[i:j],
      nrow = length(x1), ncol = length(x2),
      byrow = TRUE
    )

    se <- matrix(se[i:j],
      nrow = length(x1), ncol = length(x2),
      byrow = TRUE
    )

    perc1 <- anova_percent[x1_name, "y"]
    perc2 <- anova_percent[x2_name, "y"]
    perc12 <- anova_percent[paste(x1_name, ":", x2_name,
      sep = ""
    ), "y"]
    tot.perc <- perc1 + perc2 + perc12
    # if (is.na(tot.perc))
    # perc1 or perc2 are NA.
    # tot.perc <- perc12


    plot_text <- paste(y_name, "(", x1_name, ", ", x2_name, ") : ",
      .ANOVAFormat(perc1), "+", .ANOVAFormat(perc2),
      "+", .ANOVAFormat(perc12),
      "=", .ANOVAFormat(tot.perc), "%",
      sep = ""
    )

    resize <- .X1X2Resize(x1, x2, x1_name, x2_name, x_units, x_names)
    x1 <- resize$x1
    x2 <- resize$x2
    x1_lab <- resize$x1_lab
    x2_lab <- resize$x2_lab
    .PlotContour(x1, x2, y, y_values, x1_lab, x2_lab, plot_text)


    plot_text <- paste("se[", y_name, "(", x1_name, ", ",
      x2_name, ")]",
      sep = ""
    )


    # if (title) {
    #   title(paste(
    #     "Joint Effects for", y_name, y_units
    #   ))
    # }

    if (se_plot) {
      .PlotContour(
        x1, x2, se, se_values, x1_lab, x2_lab,
        plot_text
      )
    }

    # # Has to be repeated in case se is on a new page.
    # title(paste("Joint Effects for", y_name, y_units))

    i <- j + 1
  }
}

.XNameExtract <- function(anova_percent) {
  anova_rnames <- rownames(anova_percent)
  x_names <- NULL
  for (i in 1:length(anova_rnames)) {
    if (!grepl(":", anova_rnames[i], fixed = TRUE)) {
      x_names <- c(x_names, anova_rnames[i])
    } else {
      return(x_names)
    }
  }
}

.XUnitsFind <- function(x, x_names) {
  for (i in 1:length(x_names)) {
    if (x == x_names[i]) {
      return(i)
    }
  }
}

.X1X2Resize <- function(x1, x2, x1_name, x2_name, x_units, x_names) {
  i <- 0
  j <- 0
  while (max(abs(x1)) < 0.01) {
    x1 <- x1 * 1000.0
    i <- i + 3
  }
  while (max(abs(x1)) > 100.0) {
    x1 <- x1 / 1000.0
    i <- i - 3
  }
  while (max(abs(x2)) < 0.01) {
    x2 <- x2 * 1000.0
    j <- j + 3
  }
  while (max(abs(x2)) > 100.0) {
    x2 <- x2 / 1000.0
    j <- j - 3
  }
  x1_lab <- NULL
  x2_lab <- NULL
  if (!is.null(x_units)) {
    if (i == 0) {
      x1_lab <- paste(x1_name, " ", "(", x_units[.XUnitsFind(x1_name, x_names)], ")", sep = "") # Format: x1_name (x1_unit)
    } else if (i > 0) {
      x1_lab <- paste0(x1_name, " ", "(10^(", -i, ") ", x_units[.XUnitsFind(x1_name, x_names)], ")", sep = "") # Format: x1_name (10^(-n) x1_unit)
    } else {
      x1_lab <- paste0(x1_name, " ", "(10^", -i, " ", x_units[.XUnitsFind(x1_name, x_names)], ")", sep = "") # Format: x1_name (10^n x1_unit)
    }
    if (j == 0) {
      x2_lab <- paste(x2_name, " ", "(", x_units[.XUnitsFind(x2_name, x_names)], ")", sep = "") # Format: x1_name (x1_unit)
    } else if (j > 0) {
      x2_lab <- paste0(x2_name, " ", "(10^(", -j, ") ", x_units[.XUnitsFind(x2_name, x_names)], ")", sep = "") # Format: x1_name (10^(-n) x1_unit)
    } else {
      x2_lab <- paste0(x2_name, " ", "(10^", -j, " ", x_units[.XUnitsFind(x2_name, x_names)], ")", sep = "") # Format: x1_name (10^n x1_unit)
    }
  } else {
    x1_lab <- x1_name
    x2_lab <- x2_name
  }
  scaled_result <- NULL
  scaled_result$x1 <- x1
  scaled_result$x2 <- x2
  scaled_result$x1_lab <- x1_lab
  scaled_result$x2_lab <- x2_lab
  scaled_result
}

.PlotContour <- function(x1, x2, y, contour_values = NULL,
                         x1_lab = "x1", x2_lab = "x2", plot_text = "") {
  if (is.null(contour_values)) {
    suppressWarnings(
      contour(x1, x2, y, xlab = x1_lab, ylab = x2_lab)
    )
  } else {
    suppressWarnings(
      contour(x1, x2, y,
        levels = contour_values,
        xlab = x1_lab, ylab = x2_lab
      )
    )
  }

  title(plot_text)
}


.ANOVAFormat <- function(x)
                         # Return a character string representing x, rounded to always have
# one decimal place, possibly "0".
{
  if (is.na(x)) {
    return("NA")
  }
  x_round <- round(x, 1)
  if (x_round == round(x_round, 0)) {
    return(paste0(x_round, ".0", sep = ""))
  } else {
    return(paste0(x_round))
  }
}
