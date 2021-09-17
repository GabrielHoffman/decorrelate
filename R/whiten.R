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
#' library(Rfast)
#' set.seed(1)
#' n = 800 # number of samples
#' p = 200 # number of features
#' 
#' # create correlation matrix
#' Sigma = autocorr.mat(p, .9)
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