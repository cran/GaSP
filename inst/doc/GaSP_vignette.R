## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)

## ----setup--------------------------------------------------------------------
library(GaSP)
x <- borehole$x
y <- borehole$y
x_pred <- borehole$x_pred
y_true <- borehole$y_true

## ----data---------------------------------------------------------------------
head(borehole$x, n = 3)
head(borehole$y, n = 3)

## ----reg_model, message=FALSE, warning=FALSE, results='hide'------------------
reg_model = ~ 1

## ----reg_model_first, message=FALSE, warning=FALSE, results='hide'------------
reg_model_first = ~ 1 + r + rw + Tu

## ----reg_model_bizarre, message=FALSE, warning=FALSE, results='hide'----------
reg_model_bizarre <- ~ 1 + (r + rw + Tu)^2 + I(Hu^2)

## ----sp_model_three_inputs, message=FALSE, warning=FALSE, results='hide'------
sp_model = ~ r + rw + Tu

## ----sp_model_bizarre, message=FALSE, warning=FALSE, results='hide'-----------
sp_model_bizarre = ~ (r + rw + Tu)^2 + I(Hu^2)

## ----GaSPModel, message=FALSE, warning=FALSE, results='hide'------------------
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
borehole_gasp <- GaSPModel(
  x = x, y = y,
  reg_model = ~1, cor_family = "PowerExponential",
  cor_par = cor_par, random_error = FALSE,
  sp_var = sp_var
)

## ----Fit_mle, message=FALSE, warning=FALSE, results='hide'--------------------
borehole_fit <- Fit(
  reg_model = ~1, x = x, y = y, cor_family = "PowerExponential",
  random_error = TRUE, fit_objective = "Likelihood", model_comparison = "Objective"
)

## ----Fit_map, message=FALSE, warning=FALSE, results='hide'--------------------
borehole_fit <- Fit(
  reg_model = ~1, x = x, y = y, cor_family = "Matern",
  random_error = FALSE, nugget = 0, fit_objective = "Posterior"
)

## ----Fit_warning, error=TRUE, warning=TRUE------------------------------------
borehole_fit <- Fit(x = x, y = y, 
  reg_model = ~ 1 + a, sp_model = ~ 1 + r, random_error = FALSE
)

## ----data_pred_res------------------------------------------------------------
head(borehole$x_pred, n = 3)
head(borehole$y_true, n = 3)

## ----Pred, message=FALSE, warning=FALSE, results='hide'-----------------------
borehole_pred <- Predict(
  GaSP_model = borehole_gasp,
  x_pred = x_pred,
  generate_coefficients = TRUE
)

## ----data_pred----------------------------------------------------------------
head(borehole_pred$y_pred, n = 3)

## ----CV, message=FALSE, warning=FALSE-----------------------------------------
borehole_cv <- CrossValidate(borehole_gasp)
head(borehole_cv, n = 3)

## ----figures-p1, fig.show="hold", out.width="50%"-----------------------------
PlotPredictions(borehole_pred$y_pred, y_true,
  y_name = "Water Flow Rate", y_units = "m^3/yr", title = "Predict")
PlotStdResiduals(borehole_pred$y_pred, y_true,
  y_name = "Water Flow Rate", y_units = "m^3/yr", title = "Predict")

## ----figures-p2, fig.show="hold", out.width="50%"-----------------------------
PlotResiduals(x_pred[, 1:4], borehole_pred$y_pred,
              y_true, y_name = "Water Flow Rate", y_units = "m^3/yr")

## ----figures-p3, fig.show="hold", out.width="50%", fig.align = 'center'-------
PlotQQ(borehole_pred$y_pred, y_true, y_name = "Water Flow Rate")

## ----figures-cv, results='hide', fig.show="hide"------------------------------
PlotPredictions(borehole_cv, y,
  y_name = "Water Flow Rate", y_units = "m^3/yr", title = "CrossValidate")
PlotStdResiduals(borehole_cv, y,
  y_name = "Water Flow Rate", y_units = "m^3/yr", title = "CrossValidate")
PlotResiduals(x, borehole_cv, y, y_name = "Water Flow Rate", y_units = "m^3/yr")
PlotQQ(borehole_cv, y, y_name = "Water Flow Rate")

## ----RMSE, message=FALSE, warning=FALSE---------------------------------------
RMSE(borehole_pred$y_pred$Pred, y_true, normalized = FALSE)
RMSE(borehole_pred$y_pred$Pred, y_true, normalized = TRUE)

## ----DecribeX, message=FALSE, warning=FALSE, results='hide'-------------------
borehole_x_names <- colnames(x)
borehole_min <- c(0.05, 100.00, 63070.00, 990.00, 63.10, 700.00, 1120.00, 9855.00)
borehole_max <- c(0.15, 50000.00, 115600.00, 1110.00, 116.00, 820.00, 1680.00, 12045.00)
borehole_x_desc <- DescribeX(borehole_x_names, borehole_min, borehole_max)

## ----data_x_desc--------------------------------------------------------------
borehole_x_desc

## ----Visualize, message=FALSE, warning=FALSE, results='hide'------------------
borehole_vis <- Visualize(borehole_gasp, borehole_x_desc)

## ----anova--------------------------------------------------------------------
head(borehole_vis$anova_percent, n = 3)
tail(borehole_vis$anova_percent, n = 3)

## ----Visualize_trunc, message=FALSE, warning=FALSE----------------------------
borehole_vis <- Visualize(borehole_gasp, borehole_x_desc, 
                          main_percent = 1, interaction_percent = 1)

## ----data_vis_res-------------------------------------------------------------
head(borehole_vis$main_effect, n = 3)
head(borehole_vis$joint_effect, n = 3)

## ----figures-vis_main, fig.show="hold", out.width="33.3%"---------------------
PlotMainEffects(borehole_vis$main_effect, borehole_vis$anova_percent)

## ----figures-vis_joint, fig.show="hold", out.width="33.3%"--------------------
PlotJointEffects(borehole_vis$joint_effect, borehole_vis$anova_percent)

## ----figures-all, fig.show='hide'---------------------------------------------
PlotAll(borehole_gasp, borehole_cv, borehole_vis)

