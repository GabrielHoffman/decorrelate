#' Evaluate the log determinant
#'
#' Evaluate the log determinant of the matrix
#'
#' @param Sigma.eclairs estimate of covariance/correlation matrix from \code{eclairs()} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param alpha exponent to be applied to eigen-values
#'
#' @export
logDet <- function(Sigma.eclairs, alpha = 1) {
  # first component uses non-zero eigen-values
  ld1 <- sum(alpha * with(Sigma.eclairs, log((1 - lambda) * dSq + lambda * nu)))

  # there are p-k sample eigen-values that are zero
  ld2 <- alpha * with(Sigma.eclairs, (p - k) * log(lambda * nu))

  ld1 + ld2
}
