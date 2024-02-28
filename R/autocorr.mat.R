#' Create auto-correlation matrix
#'
#' Create auto-correlation matrix
#'
#' @param p dimension of matrix
#' @param rho autocorrelation value
#'
#' @return auto-matrix of size p with parameter rho
#'
#' @examples
#' # Create 4x4 matrix with correlation between adjacent enties is 0.9
#' autocorr.mat(4, .9)
#'
#' @export
autocorr.mat <- function(p, rho) {
  mat <- diag(p)
  rho^abs(row(mat) - col(mat))
}
