#' Data for the borehole function
#'
#' Training and test data for the borehole function; see source for background.
#'
#' @format A list with the following four data frames:
#' \describe{
#'   \item{x}{8-dimensional input for 40 training runs.}
#'   \item{y}{Output (the flow) for the 40 training runs in \code{x}.}
#'   \item{x_pred}{8-dimensional input for 1000 test runs at which to predict \code{y}.}
#'   \item{y_true}{Output for the 1000 runs in \code{x_pred}.}
#' }
#' @source \url{https://www.sfu.ca/~ssurjano/borehole.html}
"borehole"
