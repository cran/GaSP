#' Visualize a \code{GaSPModel} object.
#'
#' Carry out a functional analysis of variance (ANOVA) of a \code{GaSPModel} object
#' and generate plotting coordinates for its estimated main and 2-input joint effects.
#'
#' @inheritParams CrossValidate
#' @param x_description A data frame describing the input variables.
#' See \code{\link{DescribeX}}.
#' @param main_percent An optional minimum percentage of variation explained
#' by an input's main effect to return the effect's plotting coordinates;
#' the default of zero gives plotting coordinates for all inputs.
#' @param interaction_percent An optional minimum percentage of variation explained
#' by the interaction effect of a pair of inputs to return the plotting coordinates
#' for their joint effect (main effects plus interaction effect);
#' the default of zero gives plotting coordinates for all pairs of inputs.
#' @return A list with the following elements:
#' \item{anova_percent}{A data frame containing the ANOVA percentages
#' for all main effects and 2-input interaction effects.}
#' \item{main_effect}{A data frame with plotting coordinates for the
#' estimated main effects.}
#' \item{joint_effect}{A data frame with plotting coordinates for the
#' estimated 2-input joint effects.}
#' \item{total_percent}{Total percentage of the prediction variation
#'  accounted for by all main effects and 2-input interaction effects.}
#' \item{average}{Overall average of the prediction function,
#' averaged with respect to all inputs.}
#' \item{SE_average}{Standard error of the overall average.}
#' @details If there are many inputs, to avoid excessive plotting of
#' many trivial joint effects set \code{interaction_percent = 1} say.
#' @references Schonlau, M. and Welch, W.J. (2006),
#' "Screening the Input Variables to a Computer Model Via Analysis of Variance and Visualization",
#' in \emph{Screening: Methods for Experimentation in Industry, Drug Discovery, and Genetics,}
#' Dean. A. and Lewis, S., eds., pp. 308-327, Springer, New York,
#' doi:10.1007/0-387-28014-6_14.
#'
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
#' }
#' borehole_vis <- Visualize(borehole_fit, borehole_x_desc)
#' @export
Visualize <- function(GaSP_model, x_description,
                      main_percent = 0,
                      interaction_percent = 0) {
  backup_options <- options()
  on.exit(options(backup_options))
  if (!inherits(GaSP_model, 'GaSPModel')) {
    stop("Not a GaSP Model.")
  }
  Check <- .ErrorSummarizer()
  validate <- .GaSPModelValidate(GaSP_model, Check)
  x <- data.frame(lapply(validate$x, as.double))
  b <- .DescribeXValidate(x_description, GaSP_model$x, Check)
  .EffectPtgCheck(main_percent, interaction_percent, Check)
  .ErrorOut(Check)
  if (!is.null(x_description$NumberLevels)) {
    x_description$NumberLevels <- as.integer(x_description$NumberLevels)
  }
  if (!b) {
    x_description$NumberLevels <- NULL
  }
  y <- as.double(validate$y)
  cor_family <- .ArgtoInt(GaSP_model$cor_family)
  random_error <- validate$random_error
  reg_model <- GaSP_model$reg_model
  sp_model <- GaSP_model$sp_model
  sp_var <- GaSP_model$sp_var
  error_var <- GaSP_model$error_var
  cor_par <- validate$cor_par
  result <- .Call(
    "visualize", x, y, reg_model, sp_model, cor_family, random_error, sp_var, error_var, cor_par,
    main_percent, interaction_percent, x_description
  )
  vis <- NULL
  vis$anova_percent <- result[[1]]
  vis$main_effect <- result[[2]]
  vis$joint_effect <- result[[3]]
  summary <- result[[4]]
  vis$total_percent <- summary[1]
  vis$average <- summary[2]
  vis$SE_average <- summary[3]
  vis
}
