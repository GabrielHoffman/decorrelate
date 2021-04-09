# Gabriel Hoffman
# March 12, 2021
#
# Eclairs:
# Estimate covariance/correlation with low rank and shrinkage


#' Class eclairs
#'
#' Class \code{eclairs} 
#'
#' @details Object storing:
#' \itemize{
#'  \item{U: }{orthonormal matrix with k columns representing the low rank component}
#'  \item{dSq: }{eigen-values so that \eqn{U diag(d^2) U^T} is the low rank component}
#'  \item{lambda: }{shrinkage parameter \eqn{\lambda} for the scaled diagonal component}
#'  \item{nu: }{diagonal value, \eqn{\nu}, of target matrix in shrinkage}
#'  \item{n: }{number of samples (i.e. rows) in the original data}
#'  \item{p: }{number of features (i.e. columns) in the original data}
#'  \item{k: }{rank of low rank component}
#' }
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
          
        cat("       Estimate covariance/correlation with low rank and shrinkage\n\n")

        cat("  Original data dims:", x$n, 'x',x$p, "\n")
        cat("  Low rank component:", length(x$dSq), "\n")
      	cat("  lambda:            ", format(x$lambda, digits=3), "\n")
      	cat("  nu:                ", format(x$nu, digits=3), "\n")
})


#' Get full covariance/correlation matrix from \link{eclairs}
#'
#' Get full covariance/correlation matrix from \link{eclairs} decomposition
#'
#' @param Sigma.eclairs eclairs decomposition
#' @param ... other arguments
#'
#' @return p x p covariance/correlation matrix
#'
#' @details The full matrix is computationally expensive to compute and uses a lot of memory for large p.  So it is better to use \link{decorrelate} or \link{mult_eclairs} to perform projections in \eqn{O(np)} time.
#'
#' @rdname getCov
#' @export
setGeneric("getCov", function(Sigma.eclairs,...) standardGeneric("getCov"))

#' @rdname getCov
#' @export
setGeneric("getCor", function(Sigma.eclairs,...) standardGeneric("getCor"))



#' @rdname getCov
#' @export
setMethod('getCov', c(Sigma.eclairs = "eclairs"),
	function(Sigma.eclairs,...){

	# if disableWarn is specified and is TRUE, set disableWarn = TRUE
	# else FALSE
	lst = list(...)
	disableWarn = FALSE
	if( !is.null(lst[['disableWarn']]) ){
		disableWarn = lst[['disableWarn']]
	}

	if( !disableWarn & Sigma.eclairs$nu == 1){
		warning("eclairs estimated correlation, so the correlation matrix is returned")
	}

	# A = x$U %*% diag(x$dSq *(1-x$lambda)) %*% t(x$U) + diag(x$lambda, x$p)
	Sigma.eclairs$U %*% ((Sigma.eclairs$dSq *(1-Sigma.eclairs$lambda)) * t(Sigma.eclairs$U)) + diag(Sigma.eclairs$nu*Sigma.eclairs$lambda, Sigma.eclairs$p)
})


#' @importFrom stats cov2cor
#' @rdname getCov
#' @export
setMethod('getCor', c(Sigma.eclairs = "eclairs"), 
	function(Sigma.eclairs,...){

	cov2cor( getCov( Sigma.eclairs, disableWarn=TRUE ) )
})






#' Estimate covariance/correlation with low rank and shrinkage
#'
#' Estimate the covariance/correlation between columns as the weighted sum of a low rank matrix and a scaled identity matrix.  The weight acts to shrink the sample correlation matrix towards the identity matrix or the sample covariance matrix towards a scaled identity matrix with constant variance.  An estimate of this form is useful because it is fast, and enables fast operations downstream.
#'
#' @param X data matrix with n samples as rows and p features as columns
#' @param k the rank of the low rank component  
#' @param lambda shrinkage parameter. If not specified, it is estimated from the data.
#' @param compute compute the 'covariance' (default) or 'correlation'
# @param center center columns of X (default: TRUE) 
# @param scale scale columns of X (default: TRUE) 
#' @param warmStart result of previous SVD to initialize values
#' @param fastLambda use fast approximation to estimate lambda (default: TRUE)
#'
#' @return \link{eclairs} object storing:
#' \itemize{
#'  \item{U: }{orthonormal matrix with k columns representing the low rank component}
#'  \item{dSq: }{eigen-values so that \eqn{U diag(d^2) U^T} is the low rank component}
#'  \item{lambda: }{shrinkage parameter \eqn{\lambda} for the scaled diagonal component}
#'  \item{nu: }{diagonal value of target matrix in shrinkage}
#'  \item{n: }{number of samples (i.e. rows) in the original data}
#'  \item{p: }{number of features (i.e. columns) in the original data}
#'  \item{k: }{rank of low rank component}
#' }
#'
#' @details
#' Compute \eqn{U}, \eqn{d^2} to approximate the covariance/correlation matrix between columns of data matrix X by \eqn{U diag(d^2 (1-\lambda)) U^T + diag(\nu * \lambda)}.  When computing the covariance matrix \eqn{\nu} is the constant variance which is the mean of all feature-wise variances.  When computing the correlation matrix, \eqn{\nu = 1}.   
#'
#' @importFrom Rfast standardise colVars
# @importFrom irlba irlba
#' @importFrom PRIMME svds
#' @importFrom methods new
#'
#' @export
eclairs = function(X, k, lambda=NULL, compute=c("covariance", "correlation"), warmStart=NULL, fastLambda = TRUE){

	stopifnot(is.matrix(X))
	compute = match.arg(compute)

	# check value of lambda
	if( ! missing(lambda) & !is.null(lambda) ){
		if( lambda < 0 || lambda > 1){
			stop("lambda must be in (0,1)")
		}
	}

	n = nrow(X)
	p = ncol(X)

	center = TRUE
	# scale=TRUE
	# if correlation is computed, scale features
	scale = ifelse( compute == "correlation", TRUE, FALSE)

	# Estimate num, the scale of the target matrix
	# needs to be computed here, because X is overwritten
	nu = ifelse( compute == "correlation", 1, mean(colVars(X)))



	# scale so that cross-product gives correlation
	# features are *columns*	
	# Standarizing sets the varianes of each column to be equal
	#	so the using a single sigSq across all columns is a good
	# 	approximation
	# X = scale(X) / sqrt(n-1)
	if( center | scale ){
		X = standardise(X, center=center, scale=scale) 
	}

	# always divide by sqrt(n-1) so that var(x[,1]) 
	# is crossprod(x[,1])/(n-1) when center is TRUE
	X = X / sqrt(n-1)

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
		# Setting nv = nu = k doesn't work with warm start
		# see https://github.com/bwlewis/irlba/issues/58
		# # but setting nu=k+1 works
		# dcmp = irlba(X, nv=k, nu=k+1, v = warmStart)
		# dcmp$u = dcmp$u[,1:k]
		# dcmp$d = dcmp$d[1:k]
		# dcmp = irlba(X, nv=k, nu=k)
		if( is.null(warmStart) ){
			dcmp = svds(X, k, isreal=TRUE)
		}else{			
			dcmp = svds(X, k, v0=warmStart$U, isreal=TRUE)
		}
	}else{
		dcmp = svd(X) 
		dcmp$u = dcmp$u[,1:k, drop=FALSE]
		dcmp$v = dcmp$v[,1:k, drop=FALSE]
		dcmp$d = dcmp$d[1:k]
	}

	# Estimate lambda 
	#################

	if( missing(lambda) | is.null(lambda) ){

		if( fastLambda ){
			# fast approximation
			lambda = est_lambda_ev( dcmp$d^2, n, p)
		}else{
			# if SVD is a full rank approximation
			if( k == min(dim(X)) ){			
				lambda = shrinkcovmat.equal_lambda( t(X) )$lambda_hat
			}else{

				# if SVD is a low rank approximation
				# estimate lambda for shrinkage
				# but use lambda reconstructed from low rank decomp
				# lambda = mvIC:::shrinkcovmat.equal_lambda( t(X) )$lambda
				X_reconstruct = with(dcmp, u %*% (d * t(v)))
				lambda = shrinkcovmat.equal_lambda( t(X_reconstruct) )$lambda_hat
			}
		}

		lambda = min(1, max(1e-6, lambda))
	}

	result = list(	U 		= dcmp$v, 
					dSq 	= dcmp$d^2, 
					lambda 	= lambda,
					nu 		= nu,
					n 		= n,
					p 		= p,
					k 		= k)

	new("eclairs",	result)
}



