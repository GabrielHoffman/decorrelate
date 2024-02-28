# Gabriel Hoffman April 7, 2021 Linear time drawing from multivariate normal
# and t distributions


#' Draw from multivariate normal and t distributions
#'
#' Draw from multivariate normal and t distributions using eclairs decomposition
#'
#' @param n sample size
#' @param mu mean vector
#' @param ecl covariance matrix as an \link{eclairs} object
#' @param v degrees of freedom, defaults to Inf.  If finite, uses a multivariate t distribution
#' @param seed If you want the same to be generated again use a seed for the generator, an integer number
#'
#' @details
#' Draw from multivariate normal and t distributions using eclairs decomposition.  If the (implied) covariance matrix is \eqn{p \times p}, the standard approach is \eqn{O(p^3)}. Taking advantage of the previously computed eclairs decomposition of rank \eqn{k}, this can be done in \eqn{O(pk^2)}.
#'
#' @return matrix where rows are samples from multivariate normal or t distribution where columns have covariance specified by \code{ecl}
#'
#' @examples
#' library(Rfast)
#'
#' n <- 800 # number of samples
#' p <- 200 # number of features
#'
#' # create correlation matrix
#' Sigma <- autocorr.mat(p, .9)
#'
#' # draw data from correlation matrix Sigma
#' Y <- rmvnorm(n, rep(0, p), sigma = Sigma, seed = 1)
#'
#' # perform eclairs decomposition
#' ecl <- eclairs(Y)
#'
#' # draw from multivariate normal
#' n <- 10000
#' mu <- rep(0, ncol(Y))
#'
#' # using eclairs decomposition
#' X.draw1 <- rmvnorm_eclairs(n, mu, ecl)
#'
#' # using full covariance matrix implied by eclairs model
#' X.draw2 <- rmvnorm(n, mu, getCov(ecl))
#'
#' # assess difference betweeen covariances from two methods
#' range(cov(X.draw1) - cov(X.draw2))
#'
#' # compare covariance to the covariance matrix used to simulated the data
#' range(cov(X.draw1) - getCov(ecl))
#'
#' @importFrom stats rchisq
#' @importFrom methods is
#' @importFrom Rfast matrnorm
#' @export
rmvnorm_eclairs <- function(n, mu, ecl, v = Inf, seed = NULL) {
  stopifnot(is(ecl, "eclairs"))

  p <- length(mu)

  # simulate from standard normal
  X <- matrnorm(n, p, seed = seed)

  # Create X %*% sqrt of Sigma alpha = 1/2 induces correlation, -1/2 removes
  # correlation
  X_transform <- mult_eclairs(X, ecl$U, ecl$dSq, ecl$lambda,
    ecl$nu,
    ecl$sigma,
    alpha = 1 / 2
  )

  if (is.infinite(v)) {
    # Multivariate Normal
    X_values <- X_transform + rep(mu, rep(n, p))
  } else {
    # Multivariate t

    # create multivariate t from inverse chi-square weights, X_transform
    # and mean
    w <- sqrt(v / rchisq(n, v))
    X_values <- w * X_transform + rep(mu, rep(n, p))
  }

  # return matrix
  X_values
}
