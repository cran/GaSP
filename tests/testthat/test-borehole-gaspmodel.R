library(GaSP)

load("./borehole/const/borehole_expects.RData")

x <- borehole$x
y <- borehole$y
cor_par <- borehole_expects$corpar[c("Theta", "Alpha")]
rownames(cor_par) <- colnames(x)
sp_var <- borehole_expects$summary$StochasticProcessVariance

borehole_object <- GaSPModel(
  x = x, y = y, reg_model = ~1, cor_family = "PowerExponential", cor_par = cor_par,
  random_error = FALSE, sp_var = sp_var
)
borehole_x_pred <- borehole$x_pred

borehole_pred <- Predict(
  GaSP_model = borehole_object,
  x_pred = borehole_x_pred,
  generate_coefficients = TRUE
)

borehole_cv <- CrossValidate(borehole_object)


context("Borehole-gaspmodel")

test_that("ypred", {
  expect_equal(borehole_pred$y_pred, borehole_expects$ypred, tolerance = 1e-4)
})

test_that("pred_coefs", {
  pred_coeffs <- abs(borehole_expects$pred_coeff$Coef) + 1e-7
  borehole_pred_coeffs <- abs(borehole_pred$pred_coeffs)
  expect_equal(borehole_pred_coeffs, pred_coeffs, tolerance = 1e-4, scale = pred_coeffs)
})

test_that("cv", {
  expect_equal(borehole_cv, borehole_expects$cv, tolerance = 1e-4)
})
