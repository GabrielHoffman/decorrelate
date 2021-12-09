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
#'  \item{rownames: }{sample names from the original matrix}
#'  \item{colnames: }{features names from the original matrix}
#'  \item{method: }{method used for decomposition}
#'  \item{call: }{the function call}
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

	if( x$nu == 1){		
		rho_mle = averageCorr( x, "MLE" )
		rho_eb = averageCorr( x, "EB" )
		cat("  Avg corr (MLE):    ", format(rho_mle, digits=3), "\n")
		cat("  Avg corr (EB):     ", format(rho_eb, digits=3), "\n")
	}

	cat("  logML:             ", format(x$logML, digits=1, scientific=FALSE), "\n")
})


#' Get full covariance/correlation matrix from \link{eclairs}
#'
#' Get full covariance/correlation matrix from \link{eclairs} decomposition
#'
#' @param Sigma.eclairs eclairs decomposition
#' @param lambda shrinkage parameter for the convex combination.
# @param forCols (TRUE) compute covariance between columns  If FALSE, compute covariance between rows
#' @param ... other arguments
#'
#' @return p x p covariance/correlation matrix
#'
#' @details The full matrix is computationally expensive to compute and uses a lot of memory for large p.  So it is better to use \link{decorrelate} or \link{mult_eclairs} to perform projections in \eqn{O(np)} time.
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
#' rownames(Y) = paste0("sample_", 1:n)
#' colnames(Y) = paste0("gene_", 1:p)
#' 
#' # eclairs decomposition
#' Sigma.eclairs = eclairs(Y)
#'
#' # extract covariance implied by eclairs decomposition
#' getCov(Sigma.eclairs)[1:3, 1:3]
#'
#' @rdname getCov
#' @export
setGeneric("getCov", function(Sigma.eclairs, lambda, ...) standardGeneric("getCov"))
# , forCols=TRUE

#' @rdname getCov
#' @export
setGeneric("getCor", function(Sigma.eclairs, lambda, ...) standardGeneric("getCor"))



#' @rdname getCov
#' @export
setMethod('getCov', c(Sigma.eclairs = "eclairs"),
	function(Sigma.eclairs, lambda, ...){

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

	if( missing(lambda) ){
		lambda = Sigma.eclairs$lambda
	}


	# if( forCols ){
		# covariance between columns
		# A = x$U %*% diag(x$dSq *(1-x$lambda)) %*% t(x$U) + diag(x$lambda, x$p)
	res = Sigma.eclairs$U %*% ((Sigma.eclairs$dSq *(1-lambda)) * t(Sigma.eclairs$U)) + diag(Sigma.eclairs$nu*lambda, Sigma.eclairs$p)
	# }else{
	# 	# covariance between rows
	# 	res = Sigma.eclairs$V %*% ((Sigma.eclairs$dSq *(1-lambda)) * t(Sigma.eclairs$V)) + diag(Sigma.eclairs$nu*lambda, Sigma.eclairs$n)
	# }

	res
})


#' @importFrom stats cov2cor
#' @rdname getCov
#' @export
setMethod('getCor', c(Sigma.eclairs = "eclairs"), 
	function(Sigma.eclairs, lambda,...){

	cov2cor( getCov( Sigma.eclairs, lambda = lambda, disableWarn=TRUE ) )
})






#' Estimate covariance/correlation with low rank and shrinkage
#'
#' Estimate the covariance/correlation between columns as the weighted sum of a low rank matrix and a scaled identity matrix.  The weight acts to shrink the sample correlation matrix towards the identity matrix or the sample covariance matrix towards a scaled identity matrix with constant variance.  An estimate of this form is useful because it is fast, and enables fast operations downstream.  The method is based on the Gaussian Inverse Wishart Empirical Bayes (GIW-EB) model.
#'
#' @param X data matrix with n samples as rows and p features as columns
#' @param k the rank of the low rank component  
#' @param lambda shrinkage parameter. If not specified, it is estimated from the data.
#' @param compute compute the 'covariance' (default) or 'correlation'
# @param center center columns of X (default: TRUE) 
# @param scale scale columns of X (default: TRUE) 
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
#' @details
#' Compute \eqn{U}, \eqn{d^2} to approximate the covariance/correlation matrix between columns of data matrix X by \eqn{U diag(d^2 (1-\lambda)) U^T + diag(\nu * \lambda)}.  When computing the covariance matrix \eqn{\nu} is the constant variance which is the mean of all feature-wise variances.  When computing the correlation matrix, \eqn{\nu = 1}.   
#'
#' @examples
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
#' rownames(Y) = paste0("sample_", 1:n)
#' colnames(Y) = paste0("gene_", 1:p)
#' 
#' # eclairs decomposition
#' Sigma.eclairs = eclairs(Y, compute="covariance")
#'
#' Sigma.eclairs
#'
#' @importFrom Rfast standardise colVars eachrow
#' @importFrom PRIMME svds
#' @importFrom irlba irlba
#' @importFrom methods new
#'
#' @export
eclairs = function(X, k, lambda=NULL, compute=c("covariance", "correlation"), warmStart=NULL){

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

	# save row and columns names since X is overwritten
	rn = rownames(X)
	cn = colnames(X)

	# scale=TRUE
	# if correlation is computed, scale features
	scale = ifelse( compute == "correlation", TRUE, FALSE)

	# Estimate nu, the scale of the target matrix
	# needs to be computed here, because X is overwritten
	nu = ifelse( compute == "correlation", 1, mean(colVars(X)))

	# scale so that cross-product gives correlation
	# features are *columns*	
	# Standarizing sets the varianes of each column to be equal
	#	so the using a single sigSq across all columns is a good
	# 	approximation
	# X = scale(X) / sqrt(n-1)
	# if( center | scale ){
	# 	X = standardise(X, center=center, scale=scale) 
	# }
	# always divide by sqrt(n-1) so that var(x[,1]) 
	# is crossprod(x[,1])/(n-1) when center is TRUE
	X = standardise(X, center=TRUE, scale=scale) / sqrt(n-1)

	if( missing(k) ){
		k = min(n,p)
	}else{
		# k cannot exceed n or p
		k = min(c(k, p, n))
	}

	# SVD of X to get low rank estimate of Sigma
	if( k < min(p, n)/3){
		if( is.null(warmStart) ){
			# dcmp = svds(X, k, isreal=TRUE)
			dcmp = irlba(X, k) # should be faster thatn PRIMME::svds
		}else{			
			dcmp = svds(X, k, v0=warmStart$U, isreal=TRUE)
		}
	}else{
		dcmp = svd(X) 

		# if k < min(n,p) truncate spectrum
		if( k < length(dcmp$d)){
			dcmp$u = dcmp$u[,seq_len(k), drop=FALSE]
			dcmp$v = dcmp$v[,seq_len(k), drop=FALSE]
			dcmp$d = dcmp$d[seq_len(k)]
		}
	}

	# Estimate lambda 
	#################

	if( missing(lambda) | is.null(lambda) ){

		# Estimate lambda by empirical Bayes, using nu as scale of target
		# Since data is scaled to have var 1 (instead of n), multiply by n
		res = estimate_lambda_eb( n*dcmp$d^2, n, p, nu)
	}else{
		# compute logML for specified lambda value
		res = estimate_lambda_eb( n*dcmp$d^2, n, p, nu, lambda=lambda)
	}

	# Modify sign of dcmp$v and dcmp$u so principal components are consistant
	# This is motivated by whitening:::makePositivDiagonal()
	# but here adjust both U and V so reconstructed data is correct
	values = sign(diag(dcmp$v))
	# dcmp$v = sweep(dcmp$v, 2, values, "*")
	# dcmp$u = sweep(dcmp$u, 2, values, "*")

	# faster version
	dcmp$v = eachrow(dcmp$v, values, "*")
	dcmp$u = eachrow(dcmp$u, values, "*")

	result = list(	U 		= dcmp$v, 
				dSq 		= dcmp$d^2, 
				V 		= dcmp$u,
				lambda 	= res$lambda,
				logML 	= res$logML,
				nu 		= nu,
				n 		= n,
				p 		= p,
				k 		= k,
				rownames	= rn,
				colnames	= cn,
				method 	= "svd",
				call 		= match.call())

	new("eclairs",	result)
}


#' Plot eclairs object
#'
#' Plot eclairs object
#' @param x eclairs object
#' @param y extra argument, not used
#' @param ... additional arguments
#'
#' @importFrom graphics points legend arrows par
#' @export
setMethod("plot", "eclairs", function(x, y, ...) {

	args = match.call(expand.dots=TRUE)

	col = args$col
	if( ! is(col, "character")){
		col = "deepskyblue"
	}

	main = args$main
	if( ! is(main, "character")){
		main = "Eigen-values"
	}

	# plot observed eigen-values
	plot(x$dSq, ylab="Eigen-values", main=main, ylim=c(0, max(x$dSq)))

	# plot shrunk eigen-values
	ev_shrunk = with(x, (1-lambda)*dSq + nu*lambda)
	points(ev_shrunk, col=col, pch='+')

	# add legend
	legend("topright", legend=c("Observed", "Shrinkage estimate"), fill=c("black", col),border="white", bty="n")

	# get max x and y values
	maxy = max(x$dSq)
	maxx = length(ev_shrunk)

	# add info about dataset
	legend( 0.7*maxx, 0.8*maxy, ncol = 2,
	  legend = c(expression(lambda~':'), expression(nu~':'), expression(n~':'), expression(p~':'), 
					format(x$lambda, digits=3),
	  			format(x$nu, digits=3),
	  			format(x$n, digits=0, big.mark=','),
	  			format(x$p, digits=0, big.mark=',')
	  	),
	  bty="n",
	  text.width=0	)

	# if SVD is low rank, 
	isLowRank = length(x$dSq) < x$p

	if( isLowRank ){
		# start at last eigen-value, and end at right side of plot
		xvals = c(maxx + .2, par("usr")[2])
		yvals = rep(x$lambda * x$nu,2)

		arrows( xvals[1], yvals[1], xvals[2], yvals[2], length=.1, col=col)
	}
})








