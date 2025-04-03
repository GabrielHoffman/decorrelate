#' Evaluate the log determinant
#'
#' Evaluate the log determinant of the matrix
#'
#' @param ecl estimate of covariance/correlation matrix from \code{eclairs()} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param alpha exponent to be applied to eigen-values
#'
#' @return log determinant
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
#' rownames(Y) <- paste0("sample_", seq(n))
#' colnames(Y) <- paste0("gene_", seq(p))
#'
#' # eclairs decomposition
#' ecl <- eclairs(Y)
#'
#' logDet(ecl)
#' @export
logDet <- function(ecl, alpha = 1) {
  # first component uses non-zero eigen-values
  ld1 <- sum(alpha * with(ecl, log((1 - lambda) * dSq + lambda * nu)))

  # there are p-k sample eigen-values that are zero
  ld2 <- alpha * with(ecl, (p - k) * log(lambda * nu))

  ld3 <- 2 * alpha * sum(log(ecl$sigma))

  ld1 + ld2 + ld3
}
