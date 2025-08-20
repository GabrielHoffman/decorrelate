# Gabriel Hoffman April 7, 2021 Linear time drawing from multivariate normal
# and t distributions


#' Draw from multivariate normal and t distributions
#'
#' Draw from multivariate normal and t distributions using eclairs decomposition
#'
#' @param n sample size
#' @param mu mean vector
#' @param ecl covariance matrix as an \link{eclairs} object
#' @param nu degrees of freedom.  If finite, uses a multivariate t distribution
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
#' @importFrom stats rchisq
#' @importFrom methods is
#' @importFrom Rfast matrnorm
#' @export
#' @rdname rmvnorm
rmvnorm_eclairs <- function(n, mu, ecl, seed = NULL) {

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

  X_transform + rep(mu, rep(n, p))
}

#' @export
#' @rdname rmvnorm
rmvt_eclairs <- function(n, mu, ecl, nu, seed = NULL) {

  # Multivariate normal
  X_transform <- rmvnorm_eclairs(n, 0, ecl, seed)

  if (is.infinite(nu)) {
    # Multivariate Normal
    X_values <- X_transform
  } else {
    # Multivariate t
    # create multivariate t from inverse chi-square weights, X_transform, and mean
    w <- sqrt(v / rchisq(n, nu))
    X_values <- w * X_transform + rep(mu, rep(n, p))
  }

  # return matrix
  X_values
}



#' Density of the multivariate normal and t distributions
#'
#' Density of the multivariate normal and t distributions
#'
#' @param x numeric matrix with the data. The rows correspond to observations and the columns to variables.
#' @param mu mean vector
#' @param ecl covariance matrix as an \link{eclairs} object
#' @param nu degrees of freedom.  If finite, uses a multivariate t distribution
#'
#' @return  A numerical vector with the density values calculated at each vector (row of the matrix x).
#'
#' @seealso \code{Rfast::dmvnorm()}, \code{eclairs()}
#'
#' @examples
#' library(Rfast)
#' 
#' X <- matrnorm(10000, 10)
#' mu <- colmeans(X)
#' Sigma <- cov(X)
#' ecl <- eclairs(X, lambda=1e-8)
#' 
#' # Multivariate normal
#' dmvnorm(X[1,], mu, Sigma, log=TRUE) 
#' dmvnorm_eclairs( X[1,], mu, ecl, log=TRUE)
#'
#' # Multivariate Student t
#' dmvt(X[1,], mu, Sigma, nu=5, log=TRUE) 
#' dmvt_eclairs( X[1,], mu, ecl, nu=5,  log=TRUE)
#
#' @export
#' @importFrom Rfast eachrow
#' @rdname dmvnorm
dmvnorm_eclairs <- function(x, mu, ecl, log = FALSE) {

  # Adapted from Rfast::dmvnorm
  quat <- -0.5 * mahalanobisDistance(ecl, eachrow(x, mu, '-'))
  pow <- length(mu)/2
  if( log ){
    logcon <- pow*log(2*pi) + logDet(ecl)/2
    den <- quat - logcon
  }else{
    con <- (2 * pi)^pow * sqrt(exp(logDet(ecl)))
    den <- exp(quat)/con
  }
  den
}


#' @export
#' @rdname dmvnorm
dmvt_eclairs <- function(x, mu, ecl, nu, log = FALSE) {

  # Adapted from Rfast::dmvt
  p <- length(mu)
  den <- lgamma((nu + p)/2) - lgamma(nu/2) - 0.5 * p * log(pi * 
      nu) - 0.5 * logDet(ecl) - 0.5 * (nu + p) * log1p(mahalanobisDistance(ecl, eachrow(x, mu, '-'))/nu)
  if (log) {
      den <- den
  }
  else den <- exp(den)
  den
}



