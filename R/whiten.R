# Gabriel Hoffman
# June 1, 2021
#
# whiten() function combines eclairs and decorrelate 

#' Decorrelation projection + eclairs
#' 
#' Efficient decorrelation projection using eclairs decomposition
#'
#' @param X matrix to be transformed so *columns* are independent
#' @param k the rank of the low rank component  
#' @param lambda specify lambda and override value estimated by \code{eclairs}
#'
#' @examples
#'
#' library(Matrix)
#' library(Rfast)
#' set.seed(1)
#' n = 800 # number of samples
#' p = 8*200 # number of features
#' 
#' # Create correlation matrix with autocorrelation
#' autocorr.mat <- function(p = 100, rho = 0.9) {
#' mat <- diag(p)
#' return(rho^abs(row(mat)-col(mat)))
#' }
#' 
#' # create correlation matrix
#' Sigma = autocorr.mat(p/8, .9)
#' Sigma = bdiag(Sigma, Sigma)
#' Sigma = bdiag(Sigma, Sigma)
#' Sigma = bdiag(Sigma, Sigma)
#' 
#' # draw data from correlation matrix Sigma
#' Y = rmvnorm(n, rep(0, p), sigma=Sigma*5.1)
#' 
#' # eclairs decomposition
#' ecl = eclairs(Y)
#' 
#' # whitened Y
#' Y.transform = decorrelate(Y, ecl)
#' 
#' # Combine eclairs and decorrelate into one step
#' Y.transform2 = whiten(Y)
#' 
#' @export
whiten = function(X, k = ncol(X), lambda=NULL){
  ecl = eclairs(X, k, lambda=lambda)

  decorrelate(X, ecl)
}