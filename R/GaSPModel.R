#' Create a \code{GaSPModel} object.
#'
#' Return a template for a \code{GaSPModel} object.
#'
#' @param x A data frame containing the input (explanatory variable) training data.
#' @param y A vector or a data frame with one column
#' containing the output (response) training data.
#' @param reg_model The regression model, specified as a formula,
#' but note the left-hand side of the formula is unused; see example.
#' @param sp_model An optional stochastic process model, specified as a formula,
#' but note the left-hand side of the formula and the intercept are unused.
#' The default \code{NULL} uses all column names in \code{x}.
#' @param cor_family A character string specifying the
#' (product, anisoptropic) correlation-function family:
#' "PowerExponential" for the power-exponential family
#' or "Matern" for the Matern family.
#' @param cor_par A data frame containing the correlation parameters
#' with one row per \code{sp_model} term and two columns (see Details).
#' @param random_error A boolean for the presence or not of a
#' random (measurement, white-noise) error term.
#' @param sp_var The stochastic process variance.
#' @param error_var The random error variance, with default 0.
#' @return A \code{GaSPModel} object, which is a list with the following components:
#' \item{x}{The data frame containing the input training data.}
#' \item{y}{The training output data, now as a vector.}
#' \item{reg_model}{The regression model, now in the form of a data frame.}
#' \item{sp_model}{The stochastic process model,
#' now in the form of a data frame.}
#' \item{cor_family}{The correlation family.}
#' \item{cor_par}{The data frame containing the correlation parameters.}
#' \item{random_error}{The boolean for the presence or not of a random error term.}
#' \item{sp_var}{The stochastic process variance.}
#' \item{error_var}{The random error variance.}
#' \item{beta}{A placeholder for a data frame to hold the regression-model parameters.}
#' \item{objective}{A placeholder for the maximum fit objective.}
#' \item{cond_num}{A placeholder for the condition number.}
#' \item{CVRMSE}{A placeholder for the model's
#' cross-validated root mean squared error.}
#' @details The data frame \code{cor_par} contains one row for each term
#' in the stochastic process model.
#' There are two columns.  The first is named \code{Theta},
#' and the second is either \code{Alpha} (power-exponential)
#' or \code{Derivatives} (Matern).
#' Let \eqn{h_j} be a distance between points for term \eqn{j}
#' in the stochastic-process model.
#' For power-exponential, the contribution to the product correlation
#' from term \eqn{j} depends on
#' a distance-scale parameter \eqn{\theta_j} from the \code{Theta} column and
#' a smoothness parameter \eqn{\alpha_j} from the \code{Alpha} column;
#' the contribution is \eqn{exp(-\theta_j h_j^{2 - \alpha_j})}.
#' For example, \eqn{\alpha_j = 0} gives the squared-exponential (Gaussian) correlation.
#' The contribution to the product correlation for Matern also depends on \eqn{\theta_j},
#' and the second parameter is the number of derivatives \eqn{\delta_j = 0, 1, 2, 3}
#' from the \code{Derivatives} column.
#' The contribution is
#' \eqn{exp(-\theta_j h_j)} for \eqn{\delta_j = 0} (the exponential correlation),
#' \eqn{exp(-\theta_j h_j) (\theta_j h_j + 1)} for \eqn{\delta_j = 1},
#' \eqn{exp(-\theta_j h_j) ((\theta_j h_j)^2 / 3 + \theta_j h_j + 1)} for \eqn{\delta_j = 2}, and
#' \eqn{exp(-\theta_j h_j^2)} for \eqn{\delta_j = 3} (the squared-exponential correlation).
#' Note that \eqn{\delta_j = 3} codes for a limiting infinite number of derivatives.
#' This is not the usual parameterization of the Matern, but it is
#' consistent with power-exponential for the exponential and squared-exponential
#' special cases common to both.
#'
#' A value should be given to \code{error_var} if the model
#' has a random-error term (\code{random_error = TRUE}),
#' and a small "nugget" such as \eqn{10^{-9}} may be needed for improved
#' numerical conditioning.
#'
#' @note This function does not excecute \code{\link{Fit}} and is intended
#' for \code{\link{CrossValidate}}, \code{\link{Predict}} and \code{\link{Visualize}}
#' with models trained otherwise by the user.
#' Placeholders do not need to be specified to excecute these further functions,
#' as they are always recomputed as needed.
#'
#' @references Sacks, J., Welch, W.J., Mitchell, T.J., and Wynn, H.P. (1989)
#' "Design and Analysis of Computer Experiments", \emph{Statistical Science}, 4, pp. 409-423,
#' doi:10.1214/ss/1177012413.
#' @examples
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
#' borehole_gasp <- GaSPModel(
#'   x = borehole$x, y = borehole$y,
#'   reg_model = ~1, cor_family = "PowerExponential",
#'   cor_par = cor_par, random_error = FALSE,
#'   sp_var = sp_var
#' )
#' @export
GaSPModel <- function(x, y, reg_model, sp_model = NULL,
                      cor_family = c("PowerExponential", "Matern"),
                      cor_par,
                      random_error = c(FALSE, TRUE),
                      sp_var,
                      error_var = 0) {
  backup_options <- options()
  on.exit(options(backup_options))
  cor_family <- match.arg(cor_family)
  result <- .GaSPModelCheck(
    reg_model, sp_model, x, y, cor_family, random_error,
    cor_par, sp_var, error_var
  )
  x <- result$x
  y <- result$y
  cor_par <- result$cor_par
  regressionModel <- .RegMod(reg_model)
  obj <- NULL
  obj$x <- data.frame(lapply(x, as.double))
  obj$y <- as.double(y)
  obj$reg_model <- regressionModel

  if (is.null(sp_model)) {
    sp_model_df <- NULL
    sp_model_df <- data.frame(Terms = colnames(x), stringsAsFactors = FALSE)
    obj$sp_model <- sp_model_df
  } else {
    obj$sp_model <- .SPMod(sp_model)
  }

  obj$cor_family <- cor_family
  # expanded_design matrix
  obj$cor_par <- cor_par
  obj$cor_par[[1]] <- as.double(cor_par[[1]])
  obj$cor_par[[2]] <- as.double(cor_par[[2]])

  obj$random_error <- random_error

  # sp_var
  obj$sp_var <- sp_var
  obj$error_var <- error_var

  beta_df <- NULL
  if (length(beta_df) != 0) {
    beta_df <- data.frame(Beta = c(0), stringsAsFactors = FALSE)
    if (nrow(regressionModel) > 1) {
      for (i in 2:nrow(regressionModel)) {
        beta_df[nrow(beta_df) + 1, ] <- c(0)
      }
    }
    rownames(beta_df) <- regressionModel$Terms
  } else {
    beta_df <- data.frame(NULL)
  }
  obj$beta <- beta_df

  obj$objective <- 0
  obj$cond_num <- 0
  obj$CVRMSE <- 0

  class(obj) <- "GaSPModel"
  return(obj)
}
