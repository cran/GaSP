% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Predict.R
\name{Predict}
\alias{Predict}
\title{Predict from a \code{GaSPModel} object.}
\usage{
Predict(GaSP_model, x_pred, generate_coefficients = c(FALSE, TRUE))
}
\arguments{
\item{GaSP_model}{Object of class \code{\link{GaSPModel}}.}

\item{x_pred}{A data frame containing the values of the input variables
at which to predict the output.}

\item{generate_coefficients}{A boolean indicating whether
coefficients for further external predictions are generated.}
}
\value{
A list with the following elements:
\item{y_pred}{A data frame with two columns: the predictions
\code{Pred} and their standard errors \code{SE}.}
\item{pred_coeffs}{A vector of coefficients for further predictions;
\code{NULL} if \code{generate_coefficients} is \code{FALSE}.}
}
\description{
Predict from a \code{GaSPModel} object.
}
\note{
The vector of prediction coefficients in \code{pred_coeffs}
can be used as follows.  Let \eqn{c} denote the coefficients and let
\eqn{r} denote a vector with element \eqn{i} containing the correlation
between the output at a given new point and the output at training point \eqn{i}.
Then the prediction for the output at the new point is the dot product
of \eqn{c} and \eqn{r}.

\code{\link{RMSE}} computes the root mean squared error
of the predictions.
\code{\link{PlotPredictions}} and \code{\link{PlotResiduals}}
plot the predictions or their residuals;
\code{\link{PlotStdResiduals}} and \code{\link{PlotQQ}}
plot the standardized residuals.
}
\examples{
\dontshow{
x <- borehole$x
y <- borehole$y
theta <- c(
  5.767699e+01, 0.000000e+00, 0.000000e+00, 1.433571e-06,
  0.000000e+00, 2.366557e-06, 1.695619e-07, 2.454376e-09
)
alpha <- c(
  1.110223e-16, 0.000000e+00, 0.000000e+00, 0.000000e+00,
  0.000000e+00, 0.000000e+00, 2.494862e-03, 0.000000e+00
)
cor_par <- data.frame(Theta = theta, Alpha = alpha)
rownames(cor_par) <- colnames(borehole$x)
sp_var <- 38783.7
borehole_fit <- GaSPModel(
  x = borehole$x, y = borehole$y,
  reg_model = ~1, cor_family = "PowerExponential",
  cor_par = cor_par, random_error = FALSE,
  sp_var = sp_var
)
}

borehole_pred <- Predict(
  GaSP_model = borehole_fit,
  x_pred = borehole$x_pred,
  generate_coefficients = TRUE
)
}
