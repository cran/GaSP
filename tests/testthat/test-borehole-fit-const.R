library(GaSP)

x <- borehole$x
y <- borehole$y

borehole_fit_const <- Fit(
  x = x, y = y, reg_model = ~1, cor_family = "PowerExponential",
  random_error = FALSE, nugget = 0
)
borehole_fit_const

borehole_x_pred <- borehole$x_pred

borehole_pred <- Predict(
  GaSP_model = borehole_fit_const,
  x_pred = borehole_x_pred,
  generate_coefficients = TRUE
)

borehole_rmse_const <- RMSE(borehole_pred$y_pred$Pred, borehole$y_true)

xdescrip <- DescribeX(
  x_names = colnames(x),
  x_min = c(0.05, 100, 63070, 990, 63.10, 700, 1120, 9855),
  x_max = c(0.15, 50000, 115600, 1110, 116, 820, 1680, 12045)
)

borehole_vis <- Visualize(borehole_fit_const, xdescrip, interaction_percent = 0.5)

context("Borehole-fit-const")

test_that("Objective", {
  objective_expect <- -117.882
  expect_equal(borehole_fit_const$objective, objective_expect, tolerance = 0.1)
})

test_that("CVRMSE", {
  cvrmse_expect_const <- 1.43461
  expect_equal(borehole_fit_const$CVRMSE, cvrmse_expect_const, tolerance = 0.01, scale = cvrmse_expect_const)
})

test_that("RMSE", {
  rmse_expect_const <- 0.782414
  expect_equal(borehole_rmse_const, rmse_expect_const, tolerance = 0.01, scale = rmse_expect_const)
})

test_that("ANOVA percent", {
  anova_expect_const <- c(
    8.287164e+01, 0.0, 0.0, 4.049656e+00, 0.0, 4.263654e+00, 3.900497e+00, 1.008336e+00, 0.0, 0.0, 1.166682e+00,
    0.0, 1.228731e+00, 1.085390e+00, 2.596247e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 3.720669e-05, 3.943149e-02, 1.015297e-02, 0.0, 0.0, 0.0, 5.870910e-02, 3.898349e-03, 2.225116e-02
  )
  expect_equal(borehole_vis$anova_percent$y, anova_expect_const, tolerance = 0.01)
})
