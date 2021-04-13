# Gabriel Hoffman
# April 13, 2021
#
# Create whitening matrix explicitly 

#' Get whitening matrix
#'
#' Get whitening matrix implied by \link{eclairs} decompostion
#'
#' @param Sigma.eclairs estimate of covariance/correlation matrix from \link{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param lambda specify lambda and override value from \code{Sigma.eclairs}
#'
#' @return whitening matrix

#' @examples
#' library(Rfast)
#' 
#' n = 2000
#' p = 3
#' 
#' Y = matrnorm(n,p)*10
#' 
#' # decorrelate with implicit whitening matrix
#' # give same result as explicity whitening matrix
#' ecl <- eclairs(Y, compute="covariance", lambda=0)
#' 
#' # get explicit whitening matrix
#' W = getWhiteningMatrix( ecl, lambda=0 )
#' 
#' # apply explicit whitening matrix
#' Z1 = tcrossprod(Y, W)
#' 
#' # use implicit whitening matrix
#' Z2 = decorrelate(Y, ecl)
#' 
#' range(Z1 - Z2)
#' 
#' @export
getWhiteningMatrix = function(Sigma.eclairs, lambda){

	stopifnot(is(Sigma.eclairs, "eclairs"))

	if( ! missing(lambda) ){
		Sigma.eclairs$lambda = lambda
	}

	# dSq = lambda = nu = NULL # pass R CMD check
	v = with(Sigma.eclairs, dSq*(1-lambda) + lambda*nu)

	# this corresponds to ZCA whitening
	# omit Sigma.eclairs$U for PCA whitening
	alpha = -1/2
	Sigma.eclairs$U %*% ((v^alpha)* t(Sigma.eclairs$U))
}
