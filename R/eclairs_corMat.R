# Gabriel Hoffman
# June 1, 2021
#
# extend eclairs() to run on correlation matrix

#' Estimate covariance/correlation with low rank and shrinkage
#'
#' Estimate covariance/correlation with low rank and shrinkage from the correlation matrix
#'
#' @param C sample correlation matrix between features
#' @param n number of samples used to estimate the sample correlation matrix
#' @param k the rank of the low rank component.  Defaults to the sample size, \code{n}.  
#' @param lambda shrinkage parameter. If not specified, it is estimated from the data.
#' @param warmStart result of previous SVD to initialize values
#'
#' @return \link{eclairs} object storing:
#' \itemize{
#'  \item{U: }{orthonormal matrix with k columns representing the low rank component}
#'  \item{dSq: }{eigen-values so that \eqn{U diag(d^2) U^T} is the low rank component}
#'  \item{lambda: }{shrinkage parameter \eqn{\lambda} for the scaled diagonal component}
#'  \item{nu: }{diagonal value, \eqn{\nu}, of target matrix in shrinkage}
#'  \item{n: }{number of samples (i.e. rows) in the original data}
#'  \item{p: }{number of features (i.e. columns) in the original data}
#'  \item{k: }{rank of low rank component}
#'  \item{rownames: }{sample names from the original matrix}
#'  \item{colnames: }{features names from the original matrix}
#'  \item{method: }{method used for decomposition}
#'  \item{call: }{the function call}
#' }
#'
#' @importFrom Rfast eachrow
#' @importFrom PRIMME eigs_sym
#' @importFrom irlba partial_eigen
#' @importFrom methods new
#' @importFrom Matrix diag isSymmetric
#'
#' @export
eclairs_corMat = function(C, n, k=n, lambda=NULL, warmStart=NULL){

	# check that diagonals are all equal, up to a tolerance
	if( ! all.equal(as.numeric(diag(C)), rep(1, ncol(C))) ){
		stop("All diagonal entries of C must be 1")
	}

	# check that C is symmetric
	if( ! isSymmetric(C) ){
		stop("Matrix C is not symmetric")
	}

	if( missing(k) ){
		n = k
	}

	p = nrow(C) # number of features
	nu = mean(diag(C)) # scale of identity matrix

	if( k > p/3){
		dcmp = eigen(C)
	}else{
		# partial eigen decomposition
		if( is.null(warmStart) ){
			suppressWarnings({
				dcmp <- partial_eigen(C, k)
			})
		}else{		
			dcmp <- eigs_sym(C, k, isreal=TRUE, x0=warmStart$U)
		}
	}

	# keep only positive eigen-values
	if( any(dcmp$values <= .Machine$double.eps) ){
		keep = which(dcmp$values >.Machine$double.eps)

		dcmp$values = dcmp$values[keep]
		dcmp$vectors = dcmp$vectors[keep,,drop=FALSE]
	}

	if( missing(lambda) | is.null(lambda) ){

		# Estimate lambda by empirical Bayes, using nu as scale of target
		# Since data is scaled to have var 1 (instead of n), multiply by n
		lambda = estimate_lambda_eb( n*dcmp$values, n, p, nu)
	}

	# Modify sign of dcmp$v and dcmp$u so principal components are consistant
	# This is motivated by whitening:::makePositivDiagonal()
	# but here adjust both U and V so reconstructed data is correct
	values = sign(diag(dcmp$vectors))

	# faster version
	dcmp$vectors = eachrow(dcmp$vectors, values, "*")

	result = list(	U 		= dcmp$vectors, 
					dSq 	= dcmp$values, 
					V 		= NULL,
					lambda 	= lambda,
					nu 		= nu,
					n 		= n,
					p 		= p,
					k 		= k,
					rownames= NULL,
					colnames= colnames(C),
					method 	= "eigen",
					call 	= match.call())

	new("eclairs",	result)
}

