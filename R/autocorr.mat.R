
#' Create autocorrelation matrix
#'
#' Create autocorrelation matrix
#' 
#' @param p dimension of matrix
#' @param rho autocorrelation value
#'
#' @examples
#' # Create 10x10 matrix with correlation between adjacent enties is 0.9
#' autocorr.mat( 4, .9)
#'
#' @export
autocorr.mat <- function(p, rho) {
	mat <- diag(p)
	rho^abs(row(mat)-col(mat))
}