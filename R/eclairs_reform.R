# Gabriel Hoffman
# April 12, 2021
#
# Recompute eclairs after reforming original data, and droppping features

#' Recompute eclairs after dropping features
#'
#' Recompute eclairs after dropping features
#'
#' @param Sigma.eclairs covariance/correlation matrix as an \link{eclairs} object
#' @param k the rank of the low rank component  
#' @param drop array of variable names to drop. 
#'
#' @details Reform the dataset from the eclairs decomposition, drop features, then recompute the eclairs decomposition.  If the original SVD/eigen was truncated, then the reconstruction of the original data will be approximate.  Note that the target shrinkage matrix is the same as in \code{Sigma.eclairs}, so \eqn{\nu} is not recomputed from the retained features. 
#'
#' @return \link{eclairs} decomposition for a subset of features
#'
#' @examples
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
#' rownames(Y) = paste0("sample_", 1:n)
#' colnames(Y) = paste0("gene_", 1:p)
#' 
#' # Correlation
#' #------------
#' 
#' # eclairs decomposition
#' Sigma.eclairs = eclairs(Y, compute="correlation")
#' 
#' # features to drop
#' drop = paste0("gene_",1:100)
#' 
#' # Compute SVD on subset of eclairs decomposition
#' ecl1 = reform_decomp( Sigma.eclairs, drop=drop)
#' 
#' ecl1
#' 
#' @importFrom utils head
#' @export
reform_decomp = function(Sigma.eclairs, k = Sigma.eclairs$k, drop=NULL){

	stopifnot(is(Sigma.eclairs, "eclairs"))

	if( is.null(drop) | length(drop) == 0){
		warning("No variables are dropped.  Returning original eclairs")
		return(Sigma.eclairs)
	}

	if( !(Sigma.eclairs$method %in% c("svd", "eigen")) ){
		stop("method must be 'svd' or 'eigen'")
	}

	# check if there are entries in drop that are not in Sigma.eclairs
	problematic = which(!(drop %in% Sigma.eclairs$colnames))

	# if yes, throw warning
	if( length(problematic) > 0){
		warning(length(problematic), ' entries in drop are not in Sigma.eclairs.\nThe first few entries not found are: ', paste0("'", paste(head(drop[problematic], 3), collapse="', '"), "'"), immediate.=TRUE)
	}

	# keep features that are not in drop
	keep = which(!Sigma.eclairs$colnames %in% drop)

	V = dSq = U = NULL # pass R CMD BiocCheck

	n = Sigma.eclairs$n
	p = length(keep)

	# k cannot exceed n or p
	k = min(c(k, p, n))

	# reconstruct original dataset from SVD
	if( Sigma.eclairs$method == "svd" ){

		X_reconstruct = with(Sigma.eclairs, V %*% (sqrt(dSq) * t(U)))

		ecl = eclairs( 	X_reconstruct[, keep,drop=FALSE], 
						k = k, 
						lambda = Sigma.eclairs$lambda,
						warmStart = list(U = Sigma.eclairs$U[keep,seq_len(k),drop=FALSE]))

	}else{
		# if( Sigma.eclairs$method == "eigen" )
		# Reconstruct original dataset from eigen decomposition

		C_reconstruct = with(Sigma.eclairs, U %*% (dSq * t(U)))

		ecl = eclairs_corMat( 	C_reconstruct[keep, keep,drop=FALSE],
						n = n, 
						k = k, 
						lambda = Sigma.eclairs$lambda,
						warmStart = list(U = Sigma.eclairs$U[keep,seq_len(k),drop=FALSE]))
	}

	ecl$nu 	= Sigma.eclairs$nu
	ecl$call = match.call()

	ecl
}


# reform_decomp = function(Sigma.eclairs, k = Sigma.eclairs$k, drop=NULL){

# 	stopifnot(is(Sigma.eclairs, "eclairs"))

# 	if( is.null(drop) | length(drop) == 0){
# 		warning("No variables are dropped.  Returning original eclairs")
# 		return(Sigma.eclairs)
# 	}

# 	# check if there are entries in drop that are not in Sigma.eclairs
# 	problematic = which(!(drop %in% Sigma.eclairs$colnames))

# 	# if yes, throw warning
# 	if( length(problematic) > 0){
# 		warning(length(problematic), ' entries in drop are not in Sigma.eclairs.\nThe first few entries not found are: ', paste0("'", paste(head(drop[problematic], 3), collapse="', '"), "'"), immediate.=TRUE)
# 	}

# 	# keep features that are not in drop
# 	keep = which(!Sigma.eclairs$colnames %in% drop)

# 	V = dSq = U = NULL # pass R CMD BiocCheck

# 	n = Sigma.eclairs$n
# 	p = length(keep)

# 	# k cannot exceed n or p
# 	k = min(c(k, p, n))

# 	# reconstruct original dataset from SVD
# 	if( Sigma.eclairs$method == "svd" ){

# 		X_reconstruct = with(Sigma.eclairs, V %*% (sqrt(dSq) * t(U)))

# 		# ecl = eclairs( X_reconstruct, k, warmStart = Sigma.eclairs)
# 		# ecl$call = match.call()

# 		if( k < min(p, n)/3){

# 			# perform iterative SVD on subset of reconstructed matrix
# 			dcmp = svds(X_reconstruct[, keep,drop=FALSE], 
# 						NSvals = k, 
# 						# u0=Sigma.eclairs$V[,seq_len(k)], 
# 						v0=Sigma.eclairs$U[keep,seq_len(k),drop=FALSE], 
# 						isreal=TRUE)
# 		}else{
# 			dcmp = svd(X_reconstruct[, keep]) 
# 			dcmp$u = dcmp$u[,seq_len(k), drop=FALSE]
# 			dcmp$v = dcmp$v[,seq_len(k), drop=FALSE]
# 			dcmp$d = dcmp$d[seq_len(k)]
# 		}

# 		# Modify sign of dcmp$v and dcmp$u so principal components are consistant
# 		# This is motivated by whitening:::makePositivDiagonal()
# 		# but here adjust both U and V so reconstructed data is correct
# 		values = sign(diag(dcmp$v))
# 		dcmp$v = sweep(dcmp$v, 2, values, "*")
# 		dcmp$u = sweep(dcmp$u, 2, values, "*")
# 	}else{
# 	}

# 	# return result as eclairs decomposition
# 	result = list(	U 		= dcmp$v, 
# 					dSq 	= dcmp$d^2, 
# 					V 		= dcmp$u,
# 					lambda 	= Sigma.eclairs$lambda,
# 					nu 		= Sigma.eclairs$nu, # should this be recomputed?
# 					n 		= Sigma.eclairs$n,
# 					p 		= p,
# 					k 		= k,
# 					rownames= Sigma.eclairs$rownames,
# 					colnames= Sigma.eclairs$colnames[keep],
# 					method 	= "svd",
# 					call 	= match.call())

# 	new("eclairs",	result)
# }

