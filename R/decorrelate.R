# Gabriel Hoffman March 4, 2021 Linear time decorrelation projection using
# eclairs estimate of correlation


#' Multiply by eclairs matrix
#'
#' Multiply by \link{eclairs} matrix using special structure to achieve linear instead of cubic time complexity.
#'
#' @param X matrix to be transformed so *columns* are independent
#' @param U1 orthonormal matrix with k columns representing the low rank component
#' @param dSq1 eigen values so that \eqn{U_1 diag(d_1^2) U_1^T} is the low rank component
#' @param lambda shrinkage parameter for the convex combination.
#' @param nu diagonal value of target matrix in shrinkage
#' @param alpha exponent to be evaluated
#' @param sigma standard deviation of each feature
#' @param transpose logical, (default FALSE) indicating if X should be transposed first
#'
#' @details
#' Let \eqn{\Sigma = U_1 diag(d_1^2) U_1^T * (1-\lambda) + diag(\nu\lambda, p)}, where \eqn{\lambda} shrinkage parameter for the convex combination between a low rank matrix and the diagonal matrix with values \eqn{\nu}.
#'
#' Evaluate \eqn{X \Sigma^\alpha} using special structure of the \link{eclairs} decomposition in \eqn{O(k^2p)} when there are \eqn{k} components in the decomposition.
#'
#' @return a matrix product
#'
# @import Rcpp
#' @useDynLib decorrelate
#' @importFrom Matrix crossprod tcrossprod t
#' @export
mult_eclairs <- function(X, U1, dSq1, lambda, nu, alpha, sigma, transpose = FALSE) {
  v <- dSq1 * (1 - lambda) + lambda * nu

  isCorrMatrix <- all(sigma == 1)

  if (transpose) {
    if (nrow(X) != nrow(U1)) {
      stop(
        "Incompatable dimensions:\nData matrix X has ", nrow(X), " columns, but low rank component has ",
        nrow(U1), " rows"
      )
    }

    if (!isCorrMatrix) {
      # scale rows by standard deviation of input cols
      # X <- diag(sigma^(2*alpha)) %*% X
      # X <- (sigma^(2 * alpha)) * X
      X <- dmult(X, sigma^(2 * alpha), "left")
    }

    # decorrelate rows
    X_U1 <- crossprod(X, U1)

    # when lambda is zero, avoid computing the second part
    part1 <- 0
    if (lambda > 0) {
      part1 <- ((X - tcrossprod(U1, X_U1)) * ((lambda * nu)^alpha))
    }

    # result <- crossprod((v^alpha) * t(U1), t(X_U1)) + part1
    dM <- dmult(U1, v^alpha, "right")
    result <- tcrossprod(dM, X_U1) + part1
  } else {
    if (ncol(X) != nrow(U1)) {
      stop(
        "Incompatable dimensions:\nData matrix X has ", ncol(X), " columns, but low rank component has ",
        nrow(U1), " rows"
      )
    }

    if (!isCorrMatrix) {
      # scale columns by standard deviation of input cols
      # X <- X %*% diag(sigma^(2*alpha))
      # X <- t(t(X) * (sigma^(2 * alpha)))
      X <- dmult(X, sigma^(2 * alpha), "right")
    }

    # decorrelate columns
    X_U1 <- X %*% U1

    # when lambda is zero, avoid computing the second part
    part1 <- 0
    if (lambda > 0) {
      part1 <- (X - tcrossprod(X_U1, U1)) * ((lambda * nu)^alpha)
    }

    # result <- X_U1 %*% ((v^alpha) * t(U1)) + part1
    dM <- dmult(U1, v^alpha, "right")
    result <- tcrossprod(X_U1, dM) + part1
  }
  result
}




#' Decorrelation projection
#'
#' Efficient decorrelation projection using \link{eclairs} decomposition
#'
#' @param X matrix to be transformed so *columns* are independent
#' @param ecl estimate of covariance/correlation matrix from \link{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param lambda specify lambda and override value from \code{ecl}
#' @param transpose logical, (default FALSE) indicating if X should be transposed first
#' @param alpha default = -1/2.  Exponent of eigen-values
#'
#' @details
#' Apply a decorrelation transform using the implicit covariance approach to avoid directly evaluating the covariance matrix
#'
#' @return a matrix following the decorrelation transformation
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
#' Y <- rmvnorm(n, rep(0, p), sigma = Sigma * 5.1, seed = 1)
#'
#' # eclairs decomposition
#' ecl <- eclairs(Y)
#'
#' # whitened Y
#' Y.transform <- decorrelate(Y, ecl)
#' #
#' @export
decorrelate <- function(X, ecl, lambda, transpose = FALSE, alpha = -1 / 2) {
  stopifnot(is(ecl, "eclairs"))

  # of not specified, use lambda from ecl
  if (missing(lambda)) {
    lambda <- ecl$lambda
  }

  # check value of lambda
  if (lambda < 0 || lambda > 1) {
    stop("lambda must be in (0,1)")
  }

  toArray <- FALSE

  # handle 1D array
  if (is.numeric(X) & !is.matrix(X)) {
    toArray <- TRUE
    transpose <- FALSE
    X <- matrix(X, nrow = 1)
  }

  # alpha = -1/2 gives the decorrelating projection (i.e. whitening)
  X.decorr <- mult_eclairs(X, ecl$U, ecl$dSq, lambda,
    nu = ecl$nu,
    alpha = alpha,
    sigma = ecl$sigma,
    transpose = transpose
  )

  # if input X, is array, convert it back to array
  if (toArray) {
    X.decorr <- c(X.decorr)
  }

  X.decorr
}
