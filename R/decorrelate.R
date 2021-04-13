# Gabriel Hoffman
# March 4, 2021
#
# Linear time decorrelation projection using eclairs estimate of correlation 


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
#' @param transpose logical, (default FALSE) indicating if X should be transposed first

#' @details
#' Let \eqn{\Sigma = U_1 diag(d_1^2) U_1^T * (1-\lambda) + diag(\nu\lambda, p)}, where \eqn{\lambda} shrinkage parameter for the convex combination between a low rank matrix and the diagonal matrix with values \eqn{\nu}.
#'
#' Evaluate \eqn{X \Sigma^\alpha} using special structure of the \link{eclairs} decomposition in \eqn{O(k^2p)} when there are \eqn{k} components in the decomposition.
#'
#' @return a matrix product
#'
#' @importFrom Matrix tcrossprod
#' @export
mult_eclairs = function(X, U1, dSq1, lambda, nu, alpha, transpose=FALSE){

	v = dSq1*(1-lambda) + lambda*nu

	if( transpose ){

		if( nrow(X) != nrow(U1) ){
			stop("Incompatable dimensions:\nData matrix X has ", nrow(X), ' columns, but low rank component has ', nrow(U1), " rows")
		}

		# decorrelate rows
		X_U1 = crossprod(X, U1)
		result = crossprod((v^alpha)* t(U1), t(X_U1)) + ((X - tcrossprod(U1,X_U1) ) *((lambda*nu)^alpha))
	}else{

		if( ncol(X) != nrow(U1) ){
			stop("Incompatable dimensions:\nData matrix X has ", ncol(X), ' columns, but low rank component has ', nrow(U1), " rows")
		}

		# decorrelate columns
		X_U1 = X %*% U1
 		result = X_U1 %*% ((v^alpha) * t(U1)) + (X - tcrossprod(X_U1,U1) ) *((lambda*nu)^alpha)
	}
	result
}

# mult_eclairs = function(X, U1, dSq1, lambda, nu, alpha, transpose=FALSE){

# 	v = dSq1*(1-lambda) + lambda*nu

# 	if( transpose ){
# 		X = t(X)
# 	}

# 	if( ncol(X) != nrow(U1) ){
# 		stop("Incompatable dimensions:\nData matrix X has ", ncol(X), ' columns, but low rank component has ', nrow(U1), " rows")
# 	}

# 	X_U1 = X %*% U1
#  	result = X_U1 %*% ((v^alpha) * t(U1)) + (X - tcrossprod(X_U1,U1) ) *((lambda*nu)^alpha)

#  	if( transpose ){
#  		result = t(result)
#  	}

# 	result
# }

# I trued to make the last lime faster by avoiding transpose  
#  by using eachrow() to scale each column.
#  This is faster then sweep(), but the original code is faster
# tcrossprod(X_U1, eachrow(U1, v^alpha, oper = "*")) + (X - tcrossprod(X_U1,U1) ) *(lambda^alpha)









#' Decorrelation projection
#' 
#' Efficient decorrelation projection using \link{eclairs} decomposition
#'
#' @param X matrix to be transformed so *columns* are independent
#' @param Sigma.eclairs estimate of covariance/correlation matrix from \link{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param lambda specify lambda and override value from \code{Sigma.eclairs}
#' @param transpose logical, (default FALSE) indicating if X should be transposed first
#'
#' @details
#'' FIX NOTATION HERE: Given a vector \eqn{x ~ N(0, \Sigma)} where \eqn{\Sigma = U diag(d_1^2) U^T + diag(\sigma^2)}, transform x so that it has an identity covariance matrix.  \eqn{\Sigma^{-0.5}x} is such a projection (see Strimmer whitening).  When \eqn{\Sigma} is \eqn{p \times p}, computing this projection naively is \eqn{O(p^3)}.  Here we take advantage of the fact that \eqn{\Sigma} is the sum of a low rank decomposition, plus a scaled identity matrix
#'
#' @return a matrix following the decorrelation transformation
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
#'  mat <- diag(p)
#'  return(rho^abs(row(mat)-col(mat)))
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
#' @export
decorrelate = function(X, Sigma.eclairs, lambda, transpose = FALSE){

	stopifnot(is(Sigma.eclairs, "eclairs"))

	# of not specified, use lambda from Sigma.eclairs
	if( missing(lambda) ){
		lambda = Sigma.eclairs$lambda
	}

	# check value of lambda
	if( lambda < 0 || lambda > 1){
		stop("lambda must be in (0,1)")
	}

	toArray = FALSE

	# handle 1D array
	if( is.numeric(X) & ! is.matrix(X) ){
		toArray = TRUE
		transpose = FALSE
		X = matrix(X, nrow=1)
	}

	# alpha = -1/2 gives the decorrelating projection (i.e. whitening)
	X.decorr = mult_eclairs(X, Sigma.eclairs$U, Sigma.eclairs$dSq, lambda, nu=Sigma.eclairs$nu, alpha = -1/2, transpose=transpose)

	# if input X, is array, convert it back to array
	if( toArray ){
		X.decorr = c(X.decorr)
	}

	X.decorr
}















