#' Evaluate quadratic form
#'
#' Evaluate quadratic form
#'
#' @param Sigma.eclairs estimate of covariance/correlation matrix from \code{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param A matrix
#' @param alpha default = -1/2.  Exponent of eigen-values
#'
#' @details
#' Evaluate quadratic form \eqn{A^T \Sigma^{2*alpha} A}
#'
#' @examples
#' library(Rfast)
#' set.seed(1)
#' n <- 800 # number of samples
#' p <- 200 # number of features
#'
#' # create correlation matrix
#' Sigma <- autocorr.mat(p, .9)
#'
#' # draw data from correlation matrix Sigma
#' Y <- rmvnorm(n, rep(0, p), sigma = Sigma * 5.1)
#'
#' # eclairs decomposition
#' ecl <- eclairs(Y)
#'
#' # return scalar
#' quadForm(ecl, Y[1, ])
#'
#' # return matrix
#' quadForm(ecl, Y[1:2, ])
#'
#' @export
quadForm <- function(Sigma.eclairs, A, alpha = -1 / 2) {
  if (is.numeric(A) & !is.matrix(A)) {
    A <- matrix(A, nrow = 1)
  }

  # quadratic form
  res <- tcrossprod(decorrelate(A, Sigma.eclairs, alpha = alpha))

  # convert to scalar
  if (length(res) == 1) {
    res <- res[1]
  }

  res
}
