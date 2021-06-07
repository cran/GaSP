#' Fit a GaSP model.
#'
#' Fit (train) a GaSP model.
#'
#' @inheritParams GaSPModel
#' @param cor_par An optional data frame containing the correlation parameters
#' with one row per \code{sp_model} term and two columns set up as
#' described in \code{\link{GaSPModel}} Details;
#' only used to start the first objective maximization (see Details).
#' @param sp_var,error_var The stochastic process and error variances;
#' legal values are only used if \code{random_error = TRUE}
#' to start the first objective maximization (see Details).
#' @param nugget For numerical stability the proportion of the total variance
#' due to random error is fixed at this value (\code{random_error = FALSE}) or
#' bounded below by it (\code{random_error = TRUE}).
#' @param tries Number of optimizations of the objective from different random
#' starting points.
#' @param seed The random-number seed to generate starting points.
#' @param fit_objective The objective that \code{Fit} attempts to maximize:
#' "Likelihood" (maximum likelihood estimation)
#' or "Posterior" (Bayesian maximum a posteriori estimation).
#' @param theta_standardized_min,theta_standardized_max
#' The minimum and maximum of the standardized \eqn{\theta} parameter (see Details).
#' @param alpha_min,alpha_max The minimum and maximum
#' of the \eqn{\alpha} parameter of power-exponential.
#' @param derivatives_min,derivatives_max The minimum and maximum
#' of the \eqn{\delta} parameter of Matern.
#' @param log_obj_tol An absolute tolerance for terminating the maximization
#' of the log of the objective.
#' @param log_obj_diff The critical value for the change in the log objective
#' for informal tests during optimization of correlation parameters.
#' No testing is done with the default of 0;
#' a larger critical value such as 2 may give a more parsimonious model.
#' @param lambda_prior The rate parameter of an exponential prior
#' for each \eqn{\theta} parameter;
#' used only if \code{fit_objective = "Posterior"}.
#' @param model_comparison The criterion used to select from multiple solutions
#' when \eqn{\code{tries} > 1}: the objective function ("Objective")
#' or leave-one-out cross validation ("CV").
#' @return A \code{GaSPModel} object, which is a list with the following components:
#' \item{x}{The data frame containing the input training data.}
#' \item{y}{The training output data, now as a vector.}
#' \item{reg_model}{The regression model, now in the form of a data frame.}
#' \item{sp_model}{The stochastic process model, now in the form of a data frame.}
#' \item{cor_family}{The correlation family.}
#' \item{cor_par}{A data frame for the estimated correlation parameters.}
#' \item{random_error}{The boolean for the presence or not of a random error term.}
#' \item{sp_var}{The estimated stochastic process variance.}
#' \item{error_var}{The estimated random error variance.}
#' \item{beta}{A data frame holding the estimated regression-model parameters.}
#' \item{objective}{The maximum value found for the objective function:
#'  the log likelihood (fit_objective = "Likelihood")
#'  or the log posterior (fit_objective = "Posterior").}
#' \item{cond_num}{The condition number.}
#' \item{CVRMSE}{The leave-one-out cross-validation root mean squared error.}
#' @details
#' Fit numerically maximizes the profile objective function
#' with respect to the correlation parameters;
#' the mean and overall variance parameters are estimated in closed form
#' given the correlation parameters.
#'
#' A \code{cor_par} data frame supplied by the user is the starting point
#' for the first optimization try.
#' If \code{random_error = TRUE},
#' then \code{sp_var} / (\code{sp_var} + \code{error_var}) is another
#' correlation parameter to be optimized;
#' \code{sp_var} and \code{error_var} values supplied by the user
#' will initialize this parameter for the first try.
#'
#' Set \code{random_error = TRUE} to estimate the variance of the
#' random (measurement, white-noise) error;
#' a small \code{nugget} error variance is for numerical stability.
#'
#' For term \eqn{j} in the stochastic-process model,
#' the estimate of \eqn{\theta_j} is constrained between
#' \code{theta_standardized_min} / \eqn{r_j^2} and
#' \code{theta_standardized_max} / \eqn{r_j^2},
#' where \eqn{r_j} is the range of term \eqn{j}.
#' Note that \code{Fit} returns unscaled estimates relating to the original, unscaled inputs.
#'
#' @references Sacks, J., Welch, W.J., Mitchell, T.J., and Wynn, H.P. (1989)
#' "Design and Analysis of Computer Experiments", \emph{Statistical Science}, 4, pp. 409-423,
#' doi:10.1214/ss/1177012413.
#' @examples
#' x <- borehole$x
#' y <- borehole$y
#' borehole_fit <- Fit(
#'   reg_model = ~1, x = x, y = y, cor_family = "Matern",
#'   random_error = FALSE, nugget = 0, fit_objective = "Posterior"
#' )
#' @useDynLib GaSP, .registration = TRUE
#' @export
Fit <- function(x, y, reg_model, sp_model = NULL,
                cor_family = c("PowerExponential", "Matern"),
                cor_par = data.frame(0),
                random_error = c(FALSE, TRUE),
                sp_var = -1,
                error_var = -1,
                nugget = 1e-9,
                tries = 10,
                seed = 500,
                fit_objective = c("Likelihood", "Posterior"),
                theta_standardized_min = 0.0,
                theta_standardized_max = .Machine$double.xmax,
                alpha_min = 0.0,
                alpha_max = 1.0,
                derivatives_min = 0,
                derivatives_max = 3,
                log_obj_tol = 1e-5,
                log_obj_diff = 0,
                lambda_prior = 0.1,
                model_comparison = c("Objective", "CV")) {
  backup_options <- options()
  on.exit(options(backup_options))
  cor_family_num <- .ArgtoInt(match.arg(cor_family))
  fit_crit_num <- .ArgtoInt(match.arg(fit_objective))
  mod_comp_num <- .ArgtoInt(match.arg(model_comparison))

  check <- .ErrorSummarizer()
  xy <- .FitGaSPModelCheck(reg_model, sp_model, x, y, random_error, sp_var, error_var, check)
  x <- data.frame(lapply(xy$x, as.double))
  y <- as.double(xy$y)
  sp_model <- .SPModHelper(sp_model, x)
  cor_par <- .FitCorParCheck(cor_par, sp_model, cor_family, check)
  .ThetaCheck(theta_standardized_min, theta_standardized_max, check)
  .AlphaCheck(alpha_min, alpha_max, check)
  .DerivativesCheck(derivatives_min, derivatives_max, check)
  .NuggetCheck(nugget, check)
  .FitElseCheck(tries, seed, log_obj_tol, log_obj_diff, lambda_prior, check)
  .ErrorOut(check)

  reg_model <- .RegMod(reg_model)


  result <- .Call(
    "fit", x, y, reg_model, sp_model,
    cor_family_num, random_error,
    cor_par, sp_var, error_var, nugget,
    tries, seed, fit_crit_num,
    theta_standardized_min, theta_standardized_max, alpha_min, alpha_max, derivatives_min, derivatives_max,
    log_obj_tol, log_obj_diff, lambda_prior, mod_comp_num
  )

  cor_par <- result[[3]]
  sp_var <- result[[1]][4]
  error_var <- result[[1]][5]
  GaSP_model <- .GaSPModelNoCheck(x, y, reg_model, sp_model, cor_family, random_error, cor_par, sp_var, error_var)
  if (length(result[[2]]) != 0) {
    GaSP_model$beta$Beta <- result[[2]]
  }
  log_lik_add <- nrow(x) * (log(2 * pi) + 1) / 2
  GaSP_model$objective <- result[[1]][1] - log_lik_add
  GaSP_model$cond_num <- result[[1]][2]
  GaSP_model$CVRMSE <- result[[1]][3]
  GaSP_model
}

.GaSPModelNoCheck <- function(x, y, reg_model, sp_model,
                              cor_family = c("PowerExponential", "Matern"),
                              random_error = c(FALSE, TRUE),
                              cor_par,
                              sp_var,
                              error_var = 0) {
  obj <- NULL
  obj$x <- x
  obj$y <- y
  obj$reg_model <- reg_model
  obj$sp_model <- sp_model
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
  if (length(reg_model) != 0) {
    beta_df <- data.frame(Beta = c(0), stringsAsFactors = FALSE)
    if (nrow(reg_model) > 1) {
      for (i in 2:nrow(reg_model)) {
        beta_df[nrow(beta_df) + 1, ] <- c(0)
      }
    }
    rownames(beta_df) <- reg_model$Terms
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
