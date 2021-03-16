# Gabriel Hoffman
# March 12, 2021
#
# Eclairs:
# Estimate correlation with low rank and shrinkage


#' Class eclairs
#'
#' Class \code{eclairs} 
#'
#' @name eclairs-class
#' @rdname eclairs-class
#' @exportClass eclairs
# setClass("eclairs", representation("list"))
setClass("eclairs", contains= "list")

# elements in list
# U =	"matrix", 
# dSq = "array", 
# lambda 	= "numeric",
# n 	= "numeric",
# p 	= "numeric",
# k 	= "numeric"))

setMethod("show", 'eclairs',
      function(object){
	      print(object)
})


setMethod("print", 'eclairs',
      function(x){
          
        cat("\tEstimate correlation with low rank and shrinkage\n\n")

        cat("  Original data dims:", x$n, 'x',x$p, "\n")
        cat("  Low rank component:", length(x$dSq), "\n")
      	cat("  lambda:", format(x$lambda, digits=3), "\n")
})



#' Get full covariance/correlation matrix from eclairs
#'
#' Get full covariance/correlation matrix from eclairs decomposition
#'
#' @param cor.est eclairs decomposition
#'
#' @return p x p covariance/correlation matrix
#'
#' @details The full matrix is computationally expensive to compute and uses a lot of memory for large p.  So it is better to use \link{decorrelate} or \link{mult_eclairs} to perform operations in \eqn{O(np)} time.
#'
#' @rdname getCov
#' @export
setGeneric("getCov", function(cor.est) standardGeneric("getCov"))

#' @rdname getCov
#' @export
setGeneric("getCor", function(cor.est) standardGeneric("getCor"))



#' @export
setMethod('getCov', c(cor.est = "eclairs"),
	function(cor.est){

	# A = x$U %*% diag(x$dSq *(1-x$lambda)) %*% t(x$U) + diag(x$lambda, x$p)
	cor.est$U %*% ((cor.est$dSq *(1-cor.est$lambda)) * t(cor.est$U)) + diag(cor.est$lambda, cor.est$p)
})


#' @importFrom stats cov2cor
#' @export
setMethod('getCor', c(cor.est = "eclairs"), 
	function(cor.est){

	cov2cor( getCov( cor.est ) )
})






#' Estimate correlation with low rank and shrinkage
#'
#' Estimate the correlation between columns as the weighted sum of a low rank matrix and a scaled identity matrix.  The weight acts to shrink the sample correlation matrix towards the identity matrix. An estimate of this form is useful because it is fast, and enables fast operations downstream.
#'
#' @param X data matrix with n samples as rows and p features as columns
#' @param k the rank of the low rank component  
#' @param lambda shrinkage parameter
#' @param center center columns of X (default: TRUE) 
#' @param scale scale columns of X (default: TRUE) 
#' @param warmStart result of previous SVD to initialize values
#'
#' @return
#' \itemize{
#'  \item{"U"}{orthonormal matrix with k columns representing the low rank component}
#'  \item{"dSq"}{eigen values so that \eqn{U diag(dSq) U^T} is the low rank component}
#'  \item{"lambda"}{shrinkage parameter for the scaled diagonal component: \eqn{diag(sigSq, p)}}
#'  \item{"n"}{number of samples (i.e. rows) in the original data}
#'  \item{"p"}{number of features (i.e. columns) in the original data}
#'  \item{"k"}{rank of low rank component}
#' }
#'
#' @description
#' Compute U, dSq to approximate the correlation matrix between columns of data matrix X by \eqn{U diag(dSq * (1-\lambda)) + diag(\lambda) }.   
#'
#' @importFrom Rfast standardise colVars
#' @importFrom irlba irlba
#' @importFrom methods new
#'
#' @export
eclairs = function(X, k, lambda, center=TRUE, scale=TRUE, warmStart=NULL){

	if( ! is.matrix(X) ){
		stop("X must be a matrix")
	}

	n = nrow(X)
	p = ncol(X)

	# scale so that cross-product give correlation
	# features are *columns*	
	# Standarizing sets the varianes of each column to be equal
	#	so the using a single sigSq across all columns is a good
	# 	approximation
	# X = scale(X) / sqrt(n-1)
	if( center | scale ){
		X = standardise(X, center=center, scale=scale) 
	}
	if( scale ){
		X = X / sqrt(n-1)
	}

	# C.cor = cor(X)
	# C.cp = crossprod(X)
	# range(C.cor - C.cp)

	if( missing(k) ){
		k = min(n,p)
	}else{
		# k cannot exceed n or p
		k = min(c(k, p, n))
	}

	# SVD of X to get low rank estimate of Sigma
	if( k < min(p, n)/2){
		# Setting nv = nu = k doesnot work with warm start
		# see https://github.com/bwlewis/irlba/issues/58
		# # but setting nu=k+1 works
		# dcmp = irlba(X, nv=k, nu=k+1, v = warmStart)
		# dcmp$u = dcmp$u[,1:k]
		# dcmp$d = dcmp$d[1:k]
		dcmp = irlba(X, nv=k, nu=k)
	}else{
		dcmp = svd(X) 
		dcmp$u = dcmp$u[,1:k]
		dcmp$v = dcmp$v[,1:k]
		dcmp$d = dcmp$d[1:k]
	}

	if( missing(lambda) | is.null(lambda) ){
		# estimate lambda for shrinkage
		# but use lambda reconstructed from low rank decomp
		# lambda = mvIC:::shrinkcovmat.equal_lambda( t(X) )$lambda
		X_reconstruct = with(dcmp, u %*% diag(d) %*% t(v))
		lambda = shrinkcovmat.equal_lambda( t(X_reconstruct) )$lambda
		lambda = min(1, max(1e-4, lambda))
	}

	dSq = dcmp$d^2

	result = list(	U 		= dcmp$v, 
					dSq 	= dSq, 
					lambda 	= lambda,
					n 		= n,
					p 		= p,
					k 		= k,
					decomp	= dcmp)

	new("eclairs",	result)
}

