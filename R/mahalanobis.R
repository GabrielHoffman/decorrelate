# Jan 2, 2024


#' Mahalanobis Distance
#'
#' Mahalanobis Distance using \code{eclairs()} decomposition
#'
#' @param Sigma.eclairs estimate of covariance/correlation matrix from \code{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param X data matrix
#' @param lambda specify lambda and override value from ‘Sigma.eclairs’
#' @param center logical: should columns be centered internally
#'
#' @return array of distances
#'
#' @details
#' Evaluate quadratic form \eqn{(X-\mu)^T \Sigma^{-1} (X-\mu)} where covariance is estimated from finite sample
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
#' # Mahalanobis distance after mean centering
#' Y_center <- scale(Y, scale = FALSE)
#' mu <- colMeans(Y)
#'
#' # Standard R method
#' a <- mahalanobis(Y, mu, cov = cov(Y))
#'
#' # distance using eclairs decomposition, no shrinage
#' b <- mahalanobisDistance(ecl, Y_center, lambda = 0)
#' range(a - b)
#'
#' # with shrinkage
#' d <- mahalanobisDistance(ecl, Y_center)
#'
#' # centering internally
#' e <- mahalanobisDistance(ecl, Y, center = TRUE)
#' range(d - e)
#' #
#' @export
mahalanobisDistance <- function(Sigma.eclairs, X, lambda, center = FALSE) {
  if (is.numeric(X) & !is.matrix(X)) {
    X <- matrix(X, nrow = 1)
  }

  # center
  X_center <- scale(X, center = center, scale = FALSE)

  # transform
  res <- decorrelate(X_center, Sigma.eclairs, lambda = lambda)

  # sum of squares for each column
  rowSums(res^2)
}
