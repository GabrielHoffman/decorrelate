# Gabriel Hoffman April 13, 2021 Create whitening matrix explicitly

#' Get whitening matrix
#'
#' Get whitening matrix implied by \link{eclairs} decompostion
#'
#' @param ecl estimate of covariance/correlation matrix from \link{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param lambda specify lambda and override value from \code{ecl}
#'
#' @return whitening matrix

#' @examples
#' library(Rfast)
#'
#' n <- 2000
#' p <- 3
#'
#' Y <- matrnorm(n, p, seed = 1) * 10
#'
#' # decorrelate with implicit whitening matrix
#' # give same result as explicity whitening matrix
#' ecl <- eclairs(Y, compute = "covariance")
#'
#' # get explicit whitening matrix
#' W <- getWhiteningMatrix(ecl)
#'
#' # apply explicit whitening matrix
#' Z1 <- tcrossprod(Y, W)
#'
#' # use implicit whitening matrix
#' Z2 <- decorrelate(Y, ecl)
#'
#' range(Z1 - Z2)
#'
#' @export
getWhiteningMatrix <- function(ecl, lambda) {
  stopifnot(is(ecl, "eclairs"))

  if (!missing(lambda)) {
    ecl$lambda <- lambda
  }

  # pass R CMD check
  U <- dSq <- lambda <- nu <- NULL

  v <- with(ecl, dSq * (1 - lambda) + lambda * nu)

  alpha <- -1 / 2

  # when lambda is zero, avoid computing the second part
  part1 <- 0
  if (ecl$lambda > 0) {
    part1 <- with(ecl, (diag(1, p) - tcrossprod(U, U)) * ((lambda * nu)^alpha))
  }

  W <- with(ecl, U %*% ((v^alpha) * t(U)) + part1)

  if (any(ecl$sigma != 1)) {
    # W %*% diag(ecl$sigma^(2*alpha))
    W <- t(t(W) * (ecl$sigma^(2 * alpha)))
  }
  W
}
