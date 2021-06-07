#' Describe the input variables.
#'
#' Describe the input variables to set up integration or summation ranges
#' for \code{Visualize}.
#' @param x_names A vector of character strings containing the names
#' of the input variables.
#' @param x_min,x_max Vectors of the same length as \code{x_names}
#' containing the minima and maxima, respectively, of the input variables.
#' @param support Optional vector of character strings of the same length
#' as \code{x_names}. Valid strings for a variable are:
#' \code{"Continuous"} (continuous between the input's \code{x_min} and \code{x_max});
#' \code{"Fixed"} (the input's \code{x_min} must equal its \code{x_max});
#' and \code{"Grid"} (which requires the next argument).
#' @param num_levels An optional vector of integers for the number of levels of each input;
#' must be present if the \code{support} argument includes \code{"Grid"}.
#' An input's number of levels is 0 if it is \code{"Continuous"},
#' 1 if it is \code{"Fixed"},
#' or \eqn{> 1} if it is \code{"Grid"} to define an equally spaced grid
#' inclusive of the input's \code{x_min} and \code{x_max}.
#' @param distribution An optional vector of character strings of the same length
#' as \code{x_names} to define the weight distributions of the input variables.
#' Valid strings are \code{"Uniform"} or \code{"Normal"}
#' (ignored for \code{"Fixed"} inputs).
#' @return A data frame with the following columns:
#' \code{Variable} (containing \code{x_names}), \code{Min} (containing \code{x_min}),
#' and \code{Max} (containing \code{x_max}),
#' plus the optional columns \code{Support} (from \code{support}),
#' \code{NumberLevels} (from \code{num_levels}), and
#' \code{Distribution} (from \code{distribution}).
#' @note Does not check against \code{\link{GaSPModel}}
#' and all characters are CASE SENSITIVE.
#' @examples
#' borehole_x_names <- colnames(borehole$x)
#' borehole_min <- c(0.05, 100.00, 63070.00, 990.00, 63.10, 700.00, 1120.00, 9855.00)
#' borehole_max <- c(0.15, 50000.00, 115600.00, 1110.00, 116.00, 820.00, 1680.00, 12045.00)
#' borehole_x_desc <- DescribeX(borehole_x_names, borehole_min, borehole_max)
#' @export
DescribeX <- function(x_names, x_min, x_max, support = NULL, num_levels = NULL, distribution = NULL) {
  backup_options <- options()
  on.exit(options(backup_options))
  Check <- .ErrorSummarizer()
  b <- .DescribeXCheck(x_names, x_min, x_max, support, num_levels, distribution, Check)
  .ErrorOut(Check)
  x_description <- data.frame(Variable = x_names, Min = x_min, Max = x_max)
  if (!is.null(support)) {
    x_description$Support <- support
  }
  if (!is.null(num_levels) && b) {
    x_description$NumberLevels <- as.integer(num_levels)
  }
  if (!is.null(distribution)) {
    x_description$Distribution <- distribution
  }
  x_description
}
