.DerivativesCheck <- function(min, max, check) {
  if (!is.double(min) || !is.double(max)) {
    stop("'derivatives_min', 'derivatives_max' must be double 0, 1, 2 or 3.", call. = FALSE)
  }
  if (min < 0) {
    .Error("'derivatives_min' must be >= 0.", check)
  }
  if (min > max) {
    .Error("'derivatives_min' must be <= 'derivatives_max'.", check)
  }
  if (min != 0 && min != 1 && min != 2 && min != 3) {
    .Error("'derivatives_min' must be double 0, 1, 2 or 3.", check)
  }
  if (max != 0 && max != 1 && max != 2 && max != 3) {
    .Error("'derivatives_max' must be double 0, 1, 2 or 3.", check)
  }
}

.ThetaCheck <- function(min, max, check) {
  if (!is.double(min) || !is.double(max)) {
    stop("'theta_standardized_min', 'theta_standardized_max' must be double type.", call. = FALSE)
  }
  if (min < 0) {
    .Error("'theta_standardized_min' must be >= 0.", check)
  }
  if (min > max) {
    .Error("'theta_standardized_min' must be <= 'theta_standardized_max'.", check)
  }
}

.AlphaCheck <- function(min, max, check) {
  if (!is.double(min) || !is.double(max)) {
    stop("'alpha_min', 'alpha_max' must be double type.", call. = FALSE)
  }
  if (min < 0) {
    .Error("'alpha_min' must be >= 0.", check)
  }
  if (max >= 2) {
    .Error("'alpha_max' must be < 2.", check)
  }
  if (min > max) {
    .Error("'alpha_min' must be <= 'alpha_max'.", check)
  }
}

.NuggetCheck <- function(nugget, check) {
  if (!is.double(nugget)) {
    stop("'nugget' must be double type.", call. = FALSE)
  }
  if (nugget < 0 || nugget > 1) {
    .Error("'nugget' must be in [0, 1].", check)
  }
}

.YCheckHelper <- function(y) {
  if (is.list(y)) {
    if (ncol(y) != 1) {
      stop("list type 'y' must have 1 column.", call. = FALSE)
    }
    y <- as.vector(y[[1]])
  }
  if (is.matrix(y)) {
    if (ncol(y) != 1) {
      stop("matrix type 'y' must have 1 column.", call. = FALSE)
    }
    y <- as.vector(y[, 1])
  }
  y
}

.XYCheck <- function(x, y, check) {
  if (!is.list(x) && !is.matrix(x)) {
    stop("'x' must be list or matrix type.", call. = FALSE)
  }
  if (is.matrix(x)) {
    x <- as.data.frame(x)
  }
  y <- .YCheckHelper(y)
  if (!is.vector(y) || length(y) == 0) {
    stop("'y' must be vector type of length > 0.", call. = FALSE)
  }
  if (length(y) != nrow(x)) {
    stop("number of rows of 'x' must equal length of 'y'.", call. = FALSE)
  }
  xy <- NULL
  xy$x <- x
  xy$y <- y
  xy
}

.CorParCheck <- function(cor_par, sp_model, cor_family, check) {
  if (is.matrix(cor_par)) {
    cor_par <- as.data.frame(cor_par)
  }
  if (!is.list(cor_par)) {
    stop("'cor_par' must be list or matrix type.", call. = FALSE)
  }
  if (length(cor_par) != 2) {
    stop("'cor_par' must have 2 columns", call. = FALSE)
  }
  if (nrow(cor_par) != nrow(sp_model)) {
    stop("number of rows of 'cor_par' must equal number of terms in 'sp_model'.", call. = FALSE)
  }
  if (!all(rownames(cor_par) == sp_model$Terms)) {
    .Error("'cor_par' row names must be the same as the 'sp_model' terms.", check)
  }
  if (cor_family == "PowerExponential") {
    if (!all(colnames(cor_par) == c("Theta", "Alpha"))) {
      .Error("'cor_par' must have the column names 'Theta' and 'Alpha'.", check)
    }
  } else {
    if (!all(colnames(cor_par) == c("Theta", "Derivatives"))) {
      .Error("'cor_par' must have the column names 'Theta' and 'Derivatives'.", check)
    }
  }
  cor_par
}

.GaSPModelCheck <- function(reg_model, sp_model, x, y,
                            cor_family = c("PowerExponential", "Matern"),
                            random_error = c(FALSE, TRUE),
                            cor_par,
                            sp_var,
                            error_var) {
  check <- .ErrorSummarizer()
  if (!is.logical(random_error)) {
    stop("'random_error' must be logical.", call. = FALSE)
  }
  if (all(random_error == c(FALSE, TRUE))) {
    stop("'random_error' must be specified.", call. = FALSE)
  }
  if (!is.numeric(c(sp_var, error_var))) {
    stop("'sp_var' and/or 'error_var' must be numeric.", call. = FALSE)
  }
  xy <- .XYCheck(x, y, check)
  x <- xy$x
  y <- xy$y
  if (sp_var < 0 || error_var < 0) {
    .Error("'sp_var' and 'error_var' must be nonnegative.", check)
  }
  if (!is.null(sp_model)) {
    if (!inherits(sp_model, 'formula')) {
      stop("invalid 'sp_model'.", call. = FALSE)
    }
    sp_terms <- terms(sp_model)
    if (length(sp_terms) == 2) {
      if (grepl("1", paste(as.character(sp_terms[[2]]), sep = "", collapse = ""))) {
        warning("intercept term in 'sp_model' will not be used.", call. = FALSE)
      }
    }
    if (length(sp_terms) == 3) {
      if (grepl("1", paste(as.character(sp_terms[[3]]), sep = "", collapse = ""))) {
        warning("intercept term in 'sp_model' will not be used.", call. = FALSE)
      }
    }
    labels <- rownames(attr(sp_terms, "factors"))
    if (length(labels) > 0) {
      labels <- gsub("I(", "", labels, fixed = T)
      labels <- gsub(")", "", labels, fixed = T)
      labels <- unlist(lapply(strsplit(labels, split = "^", fixed = T), `[[`, 1))
      if (!all(labels %in% colnames(x))) {
        .Error("components of 'sp_model' terms must be column names in 'x'.", check)
      }
    }
  }
  if (!inherits(reg_model, 'formula')) {
    stop("invalid 'reg_model'.", call. = FALSE)
  }
  labels <- rownames(attr(terms(reg_model), "factors"))
  if (length(labels) > 0) {
    labels <- gsub("I(", "", labels, fixed = T)
    labels <- gsub(")", "", labels, fixed = T)
    labels <- unlist(lapply(strsplit(labels, split = "^", fixed = T), `[[`, 1))
    if (!all(labels %in% colnames(x))) {
      .Error("components of 'reg_model' terms must be column names in 'x'.", check)
    }
  }
  sp_model <- .SPModHelper(sp_model, x)
  cor_par <- .CorParCheck(cor_par, sp_model, cor_family, check)
  .ErrorOut(check)
  result <- xy
  result$sp_model <- sp_model
  result$cor_par <- cor_par
  result
}

.FitGaSPModelCheck <- function(reg_model, sp_model,
                               x, y,
                               random_error = c(FALSE, TRUE),
                               sp_var,
                               error_var,
                               check) {
  if (!is.logical(random_error)) {
    stop("'random_error' must be logical.", call. = FALSE)
  }
  if (all(random_error == c(FALSE, TRUE))) {
    stop("'random_error' must be specified.", call. = FALSE)
  }
  if (!is.numeric(c(sp_var, error_var))) {
    stop("'sp_var' and/or 'error_var' must be numeric.", call. = FALSE)
  }
  xy <- .XYCheck(x, y, check)
  x <- xy$x
  y <- xy$y
  if (sp_var < 0 || error_var < 0) {
    if (sp_var != -1 || error_var != -1) {
      .Error("'sp_var' and 'error_var' must be nonnegative.", check)
    }
  }
  if (!is.null(sp_model)) {
    if (!inherits(sp_model, 'formula')) {
      stop("invalid 'sp_model'.", call. = FALSE)
    }
    sp_terms <- terms(sp_model)
    if (length(sp_terms) == 2) {
      if (grepl("1", paste(as.character(sp_terms[[2]]), sep = "", collapse = ""))) {
        warning("intercept term in 'sp_model' will not be used.", call. = FALSE)
      }
    }
    if (length(sp_terms) == 3) {
      if (grepl("1", paste(as.character(sp_terms[[3]]), sep = "", collapse = ""))) {
        warning("intercept term in 'sp_model' will not be used.", call. = FALSE)
      }
    }
    labels <- rownames(attr(sp_terms, "factors"))
    if (length(labels) > 0) {
      if (any(grepl("*", labels, fixed = T))) {
        stop("'sp_model' terms contains unsupported components.", call. = FALSE)
      }
      labels <- gsub("I(", "", labels, fixed = T)
      labels <- gsub(")", "", labels, fixed = T)
      labels <- unlist(lapply(strsplit(labels, split = "^", fixed = T), `[[`, 1))
      if (!all(labels %in% colnames(x))) {
        .Error("components of 'sp_model' terms must be column names in 'x'.", check)
      }
    }
  }
  if (!inherits(reg_model, 'formula')) {
    stop("invalid 'reg_model'.", call. = FALSE)
  }
  labels <- rownames(attr(terms(reg_model), "factors"))
  if (length(labels) > 0) {
    if (any(grepl("*", labels, fixed = T))) {
      stop("'reg_model' terms contains unsupported components.", call. = FALSE)
    }
    labels <- gsub("I(", "", labels, fixed = T)
    labels <- gsub(")", "", labels, fixed = T)
    labels <- unlist(lapply(strsplit(labels, split = "^", fixed = T), `[[`, 1))
    if (!all(labels %in% colnames(x))) {
      .Error("components of 'reg_model' terms must be column names in 'x'.", check)
    }
  }
  xy
}

.FitCorParCheck <- function(cor_par, sp_model, cor_family, check) {
  if (ncol(cor_par) == 1 && nrow(cor_par) == 1 && cor_par[[1]][1] == 0) {
    return(cor_par)
  } else {
    cor_par <- .CorParCheck(cor_par, sp_model, cor_family, check)
    return(cor_par)
  }
}

.FitElseCheck <- function(tries, seed, log_obj_tol, log_obj_diff, lambda_prior, check) {
  if (!is.double(tries)) {
    stop("'tries' must be double type.", call. = FALSE)
  }
  if (!is.double(seed)) {
    stop("'seed' must be double type.", call. = FALSE)
  }
  if (!is.double(log_obj_tol)) {
    stop("'log_obj_tol' must be double type.", call. = FALSE)
  }
  if (!is.double(log_obj_diff)) {
    stop("'log_obj_diff' must be double type.", call. = FALSE)
  }
  if (!is.double(lambda_prior)) {
    stop("'lambda_prior' must be double type.", call. = FALSE)
  }
  if (tries <= 0) {
    .Error("'tries' must be > 0.", check)
  }
  if (lambda_prior <= 0) {
    .Error("'lambda_prior' must be > 0.", check)
  }
  if (log_obj_tol <= 0) {
    .Error("'log_obj_tol' must be > 0.", check)
  }
  if (log_obj_diff < 0) {
    .Error("'log_obj_diff' must be >= 0.", check)
  }
}

.XPredCheck <- function(x_pred, x, check) {
  if (!is.list(x_pred) && !is.matrix(x_pred)) {
    stop("'x_pred' must be list or matrix type.", call. = FALSE)
  }
  if (is.matrix(x_pred)) {
    x_pred <- as.data.frame(x_pred)
  }
  if (!all(unlist(lapply(x_pred, is.numeric)))) {
    .Error("'x_pred' must have only numeric columns.", check)
  }
  if (!all(unlist(lapply(x_pred, is.double)))) {
    warning("'x_pred' columns coerced to type double.", call. = FALSE)
    x_pred <- lapply(x_pred, as.double)
  }
  if (nrow(x_pred) == 0) {
    .Error("'x_pred' must have at least 1 row.", check)
  }
  if (ncol(x_pred) != ncol(x)) {
    .Error("'x_pred' must have the same number of columns as 'x'.", check)
  } else {
    if (!all(colnames(x_pred) == colnames(x))) {
      .Error("'x_pred' must have the same column names as 'x'.", check)
    }
  }
  x_pred
}

.YPredCheck <- function(y_pred, y_true, normalized, y_pred_se = c(0)) {
  check <- .ErrorSummarizer()
  y_pred <- .YCheckHelper(y_pred)
  y_true <- .YCheckHelper(y_true)
  if (!is.logical(normalized)) {
    stop("'normalized' must be logical.", call. = FALSE)
  }
  if (!is.vector(y_pred) || length(y_pred) == 0) {
    stop("'y_pred' must be vector type of length > 0.", call. = FALSE)
  }
  if (!is.vector(y_true) || length(y_true) == 0) {
    stop("'y_true' must be vector type of length > 0.", call. = FALSE)
  }
  if (!is.vector(y_pred_se) || length(y_pred_se) == 0) {
    stop("'y_pred_se' must be vector type of length > 0.", call. = FALSE)
  } else {
    if (y_pred_se != c(0)) {
      y_pred_se <- .YCheckHelper(y_pred_se)
      if (length(y_pred) != length(y_pred_se)) {
        .Error("'y_pred' and 'y_pred_se' must have equal lengths.", check)
      }
    }
  }
  if (length(y_true) != length(y_pred)) {
    .Error("'y_pred' and 'y_true' must have equal lengths.", check)
  }
  .ErrorOut(check)
  v <- NULL
  v$y_pred <- y_pred
  v$y_true <- y_true
  v$y_pred_se <- y_pred_se
  v
}

.GaSPModelValidate <- function(GaSP_Model, check) {
  if (length(GaSP_Model) != 13) {
    warning("'GaSP_Model' must have 13 components.", call. = FALSE)
  }
  x <- GaSP_Model$x
  y <- GaSP_Model$y
  reg_model <- GaSP_Model$reg_model
  sp_model <- GaSP_Model$sp_model
  cor_family <- GaSP_Model$cor_family
  cor_par <- GaSP_Model$cor_par
  random_error <- GaSP_Model$random_error
  sp_var <- GaSP_Model$sp_var
  error_var <- GaSP_Model$error_var
  if (is.null(x) || is.null(y) || is.null(reg_model) || is.null(sp_model) ||
    is.null(cor_family) || is.null(cor_par) ||
    is.null(random_error) || is.null(sp_var) || is.null(error_var)) {
    .Error("'GaSP_Model' must have non-null components to run 'Predict' and 'CrossValidate'.", check)
  }
  .RegSPValidator(reg_model, sp_model, x, check)
  xy <- .FitGaSPModelCheck(~1, NULL, x, y, random_error, sp_var, error_var, check)
  if (sp_var < 0 || error_var < 0) {
    .Error("'GaSP_Model' must have nonnegative sp_var and error_var", check)
  }
  if (random_error && error_var == 0) {
    warning("'GaSP_Model' has 'random_error' == TRUE and 'error_var' == 0.", call. = FALSE)
  }
  if (!random_error && error_var > 0) {
    warning("'GaSP_Model' 'random_error' assumed to be TRUE because 'error_var' > 0.", call. = FALSE)
    random_error <- TRUE
  }

  cor_par <- .CorParCheck(cor_par, sp_model, cor_family, check)
  xy$cor_par <- cor_par
  xy$random_error <- random_error
  xy
}

.PlotCheck <- function(y, y_pred, y_name, y_units = "", x = NULL, x_units = NULL) {
  check <- .ErrorSummarizer()
  if (!is.character(y_name) || length(y_name) != 1) {
    stop("'y_name' must be 1 character string.", call. = FALSE)
  }
  if (!is.character(y_units) || length(y_units) != 1) {
    stop("'y_units' must be 1 character string.", call. = FALSE)
  }
  y <- .YCheckHelper(y)
  if (!is.list(y_pred) && !is.matrix(y_pred)) {
    stop("'y_pred' must be list or matrix type.", call. = FALSE)
  }
  if (is.matrix(y_pred)) {
    y_pred <- as.data.frame(y_pred)
  }
  if (ncol(y_pred) != 2) {
    stop("'y_pred' must have 2 columns.", call. = FALSE)
  }
  if (!all(colnames(y_pred) == c("Pred", "SE"))) {
    warning("'y_pred' must have the column names 'Pred' and 'SE'.", call. = FALSE)
  }
  if (!is.vector(y) || length(y) == 0) {
    stop("'y' must be vector type of length > 0.", call. = FALSE)
  }
  if (length(y) != length(y_pred[[1]])) {
    .Error("'y' and 'y_pred' must have equal lengths.", check)
  }
  if (!is.null(x)) {
    if (is.matrix(x)) {
      x <- as.data.frame(x)
    }
    if (!is.list(x)) {
      stop("'x' must be list or matrix type.", call. = FALSE)
    }
    if (nrow(x) != length(y)) {
      .Error("number of rows of 'x' must equal length of 'y'.", check)
    }
  }
  if (!is.null(x_units)) {
    if (is.null(x)) {
      stop("'x_units' must be NULL when 'x' is NULL.", call. = FALSE)
    }
    if (is.list(x_units)) {
      if (ncol(x_units) != 1) {
        stop("list type 'x_units' must have 1 column.", call. = FALSE)
      }
      x_units <- as.vector(x_units[[1]])
    }
    if (is.matrix(x_units)) {
      if (ncol(x_units) != 1) {
        stop("matrix type 'x_units' must have 1 column.", call. = FALSE)
      }
      x_units <- as.vector(x_units[, 1])
    }
    if (!is.vector(x_units)) {
      stop("'x_units' must be vector type.", call. = FALSE)
    }
    if (length(x_units) != ncol(x)) {
      .Error("length of 'x_units' must equal number of columns of 'x'.", check)
    }
  }
  .ErrorOut(check)
  results <- NULL
  results$y <- y
  results$y_pred <- y_pred
  results$y_name <- y_name
  results$y_units <- y_units
  results$x <- x
  results$x_units <- x_units
  return(results)
}

.RegSPValidator <- function(reg_model, sp_model, x, check) {
  if (!is.list(reg_model) || !is.list(sp_model)) {
    stop("'reg_model' and 'sp_model' must be list type.", call. = FALSE)
  }
  labels <- sp_model$Terms
  if (length(labels) > 0) {
    labels <- unlist(strsplit(labels, ":", fixed = T))
    labels <- unlist(lapply(strsplit(labels, split = "^", fixed = T), `[[`, 1))
    if (!all(labels %in% colnames(x))) {
      .Error("components of 'sp_model' terms must be column names in 'x'.", check)
    }
  }
  labels <- reg_model$Terms[-1]
  if (length(labels) > 0) {
    labels <- unlist(strsplit(labels, ":", fixed = T))
    labels <- unlist(lapply(strsplit(labels, split = "^", fixed = T), `[[`, 1))
    if (!all(labels %in% colnames(x))) {
      .Error("components of 'reg_model' terms must be column names in 'x'.", check)
    }
  }
}

.DescribeXCheck <- function(x_names, min, max, support, num_levels, distribution, check) {
  use_num_levels <- TRUE
  if (!is.vector(x_names) || !is.character(x_names) || length(x_names) == 0) {
    stop("'x_names' must be a character vector of length > 0.", call. = FALSE)
  }
  if (!is.vector(min) || length(min) == 0) {
    stop("'min' must be a double vector type of length > 0.", call. = FALSE)
  }
  if (!is.vector(max) || length(max) == 0) {
    stop("'max' must be a double vector type of length > 0.", call. = FALSE)
  }
  if (!is.double(min)) {
    stop("'min' must be a double vector type of length > 0.", call. = FALSE)
  }
  if (!is.double(max)) {
    stop("'max' must be a double vector type of length > 0.", call. = FALSE)
  }
  x_names_l <- length(x_names)
  min_l <- length(min)
  which_sum_expect <- (min_l * (min_l + 1)) / 2
  max_l <- length(max)
  if (x_names_l != min_l || x_names_l != max_l || max_l != min_l) {
    stop("Arguments must be vectors of equal lengths.", call. = FALSE)
  }
  if (!all(min <= max)) {
    stop("'min' must be <= 'max'.", call. = FALSE)
  }
  if (!is.null(support)) {
    if (is.null(num_levels) && "Grid" %in% support) {
      stop("'num_levels' must be present if there is 'Grid' in the support vector.", call. = FALSE)
    }
    if (!is.vector(support) || length(support) == 0) {
      stop("'support' must be a character vector type of > 0 length.", call. = FALSE)
    }
    if (!is.character(support)) {
      stop("'support' must be a character vector of length > 0.", call. = FALSE)
    }
    if (x_names_l != length(support)) {
      stop("Arguments must be vectors of equal length.", call. = FALSE)
    }
    grid_which <- which("Grid" == support)
    fixed_which <- which("Fixed" == support)
    cont_which <- which("Continuous" == support)
    which_sum <- sum(c(grid_which, fixed_which, cont_which))
    if (which_sum_expect != which_sum) {
      stop("'support' contains unrecognized string(s)", call. = FALSE)
    }
    if (!is.null(num_levels)) {
      if (!is.vector(num_levels) || length(num_levels) == 0) {
        stop("'num_levels' must be vector type of length > 0.", call. = FALSE)
      }
      if (!is.double(num_levels) && !is.integer(num_levels)) {
        stop("'num_levels' must be a integer vector type of length > 0.", call. = FALSE)
      }
      if (x_names_l != length(num_levels)) {
        stop("Arguments must be vectors of equal length.", call. = FALSE)
      }
      if (sum(fixed_which) > 0) {
        fixed_num_levels <- num_levels[fixed_which]
        if (!all(fixed_num_levels == 1)) {
          .Error("'num_levels' for 'Fixed' must be 1.", check)
        }
      }
      if (sum(cont_which) > 0) {
        cont_num_levels <- num_levels[cont_which]
        if (!all(cont_num_levels == 0)) {
          .Error("'num_levels' for 'Continuous' must be 0.", check)
        }
      }
      if (sum(grid_which) > 0) {
        grid_num_levels <- num_levels[grid_which]
        if (!all(grid_num_levels > 1)) {
          .Error("'num_levels' for 'Grid' must be > 1.", check)
        }
      }
    }
    if (sum(fixed_which) > 0) {
      fixed_min <- min[fixed_which]
      fixed_max <- max[fixed_which]
      if (!all(fixed_min == fixed_max)) {
        .Error("'min' must be equal to 'max' for 'Fixed'.", check)
      }
    }
  } else if (!is.null(num_levels)) {
    warning("'num_levels will be ignored as 'support' is null.", call. = FALSE)
    use_num_levels <- FALSE
  }
  if (!is.null(distribution)) {
    if (!is.vector(distribution) || length(distribution) == 0) {
      stop("'distribution' must be a vector of character strings of length > 0.", call. = FALSE)
    }
    if (!is.character(distribution)) {
      stop("'distribution' must be a vector of character strings of length > 0.", call. = FALSE)
    }
    if (x_names_l != length(distribution)) {
      stop("Arguments must be vectors of equal length.", call. = FALSE)
    }
    norm_which <- which("Normal" == distribution)
    unif_which <- which("Uniform" == distribution)

    which_sum <- sum(c(norm_which, unif_which))
    if (which_sum_expect != which_sum) {
      pos <- which_sum_expect - which_sum
      stop("'distribution' contains unrecognized string(s)", call. = FALSE)
    }
  }
  use_num_levels
}

.DescribeXValidate <- function(x_describe, x, check) {
  b <- .DescribeXCheck(x_describe$Variable, x_describe$Min, x_describe$Max, x_describe$Support, x_describe$NumberLevels, x_describe$Distribution, check)
  if (!all(x_describe$Variable == colnames(x))) {
    stop("'x_describe' Variable names must be the same as the 'x' column names in 'GaSP_Model'.", call. = FALSE)
  }
  # if (!all(unlist(lapply(lapply(x, range), `[[`, 1)) >= x_describe$Min)) {
  #   .Error("'min' in 'x_describe' must be smaller than 'x' ranges.", check)
  # }
  # if (!all(unlist(lapply(lapply(x, range), `[[`, 2)) <= x_describe$Max)) {
  #   .Error("'max' in 'x_describe' must be greater than 'x' ranges.", check)
  # }
  b
}

.EffectPtgCheck <- function(main_effect_pct, interaction_effect_pct, check) {
  if (!is.double(main_effect_pct) || !is.double(interaction_effect_pct)) {
    stop("'main_effect_pct' and 'interaction_effect_pct' must be type double.", call. = FALSE)
  }
  if (main_effect_pct < 0 || interaction_effect_pct < 0) {
    .Error("'main_effect_pct' and 'interaction_effect_pct' must be >= 0.", check)
  }
  if (main_effect_pct > 100 || interaction_effect_pct > 100) {
    .Error("'main_effect_pct' and 'interaction_effect_pct' must be <= 100.", check)
  }
}

.EffectPlotCheck <- function(anova_percent, main_effect = NULL, joint_effect = NULL,
                             y_name, y_units, x_units = NULL, x_names, check) {
  if (!is.character(y_name) || length(y_name) != 1) {
    stop("'y_name' must be 1 character string.", call. = FALSE)
  }
  if (!is.character(y_units) || length(y_units) != 1) {
    stop("'y_units' must be 1 character string.", call. = FALSE)
  }
  if (!is.list(anova_percent) && !is.matrix(anova_percent)) {
    stop("'anova_percent' must be list or matrix type.", call. = FALSE)
  }
  if (is.matrix(anova_percent)) {
    anova_percent <- as.data.frame(anova_percent)
  }
  if (!is.character(rownames(anova_percent))) {
    stop("'anova_percent' must have character row names.", call. = FALSE)
  }
  if (!all(colnames(anova_percent) == c("y"))) {
    .Error("'anova_percent' must have column name 'y'.", check)
  }

  if (!is.null(main_effect)) {
    if (!is.list(main_effect) && !is.matrix(main_effect)) {
      stop("'main_effect' must be list or matrix type.", call. = FALSE)
    }
    if (is.matrix(main_effect)) {
      main_effect <- as.data.frame(main_effect)
    }
    if (!all(colnames(main_effect) == c("Variable.x_i", "x_i", "y", "SE"))) {
      .Error("'main_effect' must have column names 'Variable.x_i', 'x_i', 'y', and 'SE'.", check)
    }
  }
  if (!is.null(joint_effect)) {
    if (!is.list(joint_effect) && !is.matrix(joint_effect)) {
      stop("'joint_effect' must be list or matrix type.", call. = FALSE)
    }
    if (is.matrix(joint_effect)) {
      joint_effect <- as.data.frame(joint_effect)
    }
    if (!all(colnames(joint_effect) == c("Variable.x_i", "Variable.x_j", "x_i", "x_j", "y", "SE"))) {
      .Error("'joint_effect' must have column names 'Variable.x_i', 'Variable.x_j', 'x_i', 'x_j', 'y', and 'SE'.", check)
    }
  }
  if (!is.null(x_units)) {
    if (is.list(x_units)) {
      if (ncol(x_units) != 1) {
        stop("list type 'x_units' must have 1 column.", call. = FALSE)
      }
      x_units <- as.vector(x_units[[1]])
    }
    if (is.matrix(x_units)) {
      if (ncol(x_units) != 1) {
        stop("matrix type 'x_units' must have 1 column.", call. = FALSE)
      }
      x_units <- as.vector(x_units[, 1])
    }
    if (!is.vector(x_units) || !is.character(x_units)) {
      stop("'x_units' must be a vector of character strings.", call. = FALSE)
    }
    if (length(x_units) != length(x_names)) {
      .Error("'x_units' must have length equal to the number of columns in 'x'.", check)
    }
  }
  result <- NULL
  result$anova_percent <- anova_percent
  result$main_effect <- main_effect
  result$joint_effect <- joint_effect
  result$x_units <- x_units
  result
}

.JointEffectAdditCheck <- function(se_plot, y_values, se_values) {
  if (!is.logical(se_plot)) {
    stop("'se_plot' must logical.", call. = FALSE)
  }
  if (!is.null(y_values)) {
    if (is.list(y_values)) {
      if (ncol(y_values) != 1) {
        stop("list type 'y_values' must have 1 column.", call. = FALSE)
      }
      y_values <- as.vector(y_values[[1]])
    }
    if (is.matrix(y_values)) {
      if (ncol(y_values) != 1) {
        stop("matrix type 'y_values' must have 1 column.", call. = FALSE)
      }
      y_values <- as.vector(y_values[, 1])
    }
    if (!is.vector(y_values) || !is.numeric(y_values)) {
      stop("'y_values' must be numeric vector type.", call. = FALSE)
    }
  }
  if (!is.null(se_values)) {
    if (is.list(se_values)) {
      if (ncol(se_values) != 1) {
        stop("list type 'se_values' must have 1 column.", call. = FALSE)
      }
      se_values <- as.vector(se_values[[1]])
    }
    if (is.matrix(se_values)) {
      if (ncol(se_values) != 1) {
        stop("matrix type 'se_values' must have 1 column.", call. = FALSE)
      }
      se_values <- as.vector(se_values[, 1])
    }
    if (!is.vector(se_values) || !is.numeric(se_values)) {
      stop("'se_values' must be numeric vector type.", call. = FALSE)
    }
  }
  result <- NULL
  result$se_plot <- se_plot
  result$y_values <- y_values
  result$se_values <- se_values
  result
}
