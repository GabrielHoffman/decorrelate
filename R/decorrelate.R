# Gabriel Hoffman
# March 4, 2021
#
# Linear time decorrelation projection using eclairs estimate of correlation 


#' Multiply by eclairs matrix
#'
#' Multiply by eclairs matrix using special structure to achieve linear instead of cubic time complexity.
#'
#' @param X matrix to be transformed so *columns* are independent
#' @param U1 orthonormal matrix with k columns representing the low rank component
#' @param dSq1 eigen values so that U1 diag(dSq1) U1^T is the low rank component
#' @param lambda shrinkage parameter for the convex combination.
#' @param alpha exponent to be evaluated

#' @details
#' Let \eqn{\Sigma = U1 diag(dSq1) U1^T * (1-\lambda) + diag(\lambda, p)}, where \eqn{\lambda} shrinkage parameter for the convex combination between a low rank matrix and the identity matrix.
#'
#' Evaluate \eqn{X \Sigma^\alpha} in linear time using the special structure of the matrix
#'
#' @export
mult_eclairs = function(X, U1, dSq1, lambda, alpha){

	v = dSq1*(1-lambda) + lambda
	X_U1 = X %*% U1
 	X_U1 %*% ((v^alpha) * t(U1))  + (X - tcrossprod(X_U1,U1) ) *(lambda^alpha)
}



#' Decorrelation projection
#' 
#' Efficient decorrelation projection using eclairs decomposition
#'
#' @param X matrix to be transformed so *columns* are independent
#' @param cor.est estimate of correlation matrix from \link{eclairs} storing \code{U}, \code{dSq}, and \code{lambda}
#'
#' @details
#'' FIX NOTATION HERE: Given a vector \eqn{x ~ N(0, \Sigma)} where \eqn{\Sigma = U diag(dSq1) U^T + diag(\sigma^2)}, transform x so that it has an identity covariance matrix.  \eqn{\Sigma^{-0.5}x} is such a projection (see Strimmer whitening).  When \eqn{\Sigma} is \eqn{p \times p}, computing this project naively is \eqn{O(p^3)}.  Here we take advantage of the fact that \eqn{\Sigma} is the sum of a low rank decomposition, plus a scaled identity matrix
#'
#' @export
decorrelate = function(X, cor.est){

	transpose = FALSE
	if( min(dim(X)) == 1 ){
		transpose = TRUE
		X = matrix(X, nrow=1)
	}

	# alpha = -1/2 gives the decorrelating projection (i.e. whitening)
	X.decorr = mult_eclairs(X, cor.est$U, cor.est$dSq, cor.est$lambda, alpha = -1/2)

	if( transpose ) X.decorr = t(X.decorr)

	X.decorr
}















