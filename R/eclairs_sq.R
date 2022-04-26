# Gabriel Hoffman
# April 25, 2022
#
# eclairs_sq: exclairs to estimate 2*cor(Y^2)





#' Compute eclairs decomp of squared correlation matrix
#'
#' Given the eclairs decomp of C, compute the eclairs decomp of C^2 
#'
#' @param Sigma.eclairs estimate of correlation matrix from \code{eclairs()} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}

#' @param rank1 use the first 'rank' singular vectors from the SVD.  Using increasing rank1 will increase the accuracy of the estimation.  But note that the computationaly complexity is O(P choose(rank, 2)), where P is the number of features in the dataset
#' @param rank2 rank of \code{eclairs()} decomposition returned
#' @param varianceFraction fraction of variance to retain after truncating trailing eigen values
#'
#' @details Consider a data matrix X_{N x P} of P features and N samples where N << P. Let the columns of X be scaled so that C_{P x P} = XX^T.  C is often too big to compute directly since it is O(P^2) and O(P^3) to invert.  But we can compute the SVD of X in O(PN^2).
#' The goal is to compute the SVD of the matrix C^2, given only the SVD of C in less than O(P^2 time).  Here we compute this SVD of C^2 in O(PN^4) time, which is tractible for small N.
#' Moreover, if we use an SVD X = UDV^T with of rank R, we can approximate the SVD of C^2 in O(PR^4) using only D and V  
#' In practice, this can be reduced to O(P (choose(R,2) + R)^2)
#'
#' @examples
#' # Compute correlations directly and using eclairs decomp
#'
#' set.seed(1)
#' n = 600 # number of samples
#' p = 100 # number of features
#' 
#' # create correlation matrix
#' Sigma = autocorr.mat(p, .9)
#' 
#' # draw data from correlation matrix Sigma
#' Y = Rfast::rmvnorm(n, rep(0, p), sigma=Sigma)
#' rownames(Y) = paste0("sample_", 1:n)
#' colnames(Y) = paste0("gene_", 1:p)
#' 
#' # correlation computed directly
#' C = cor(Y)
#' 
#' # correlation from eclairs decomposition
#' ecl = eclairs(Y, compute="cor")
#' C.eclairs = getCor(ecl, lambda=0)
#' 
#' all.equal(C, C.eclairs)
#' 
#' # Correlation of Y^2
#' #-------------------
#' 
#' # exact quadratic way
#' C = 2*cor(Y)^2 
#' 
#' # faster low rank
#' ecl2 = eclairs_sq(ecl)
#' C.eclairs = 2*getCor(ecl2, lambda=0)
#' 
#' all.equal(C.eclairs, C)
#' 
#' @return compute the eclairs of C^2
#' @importFrom Rfast eachrow colVars
#' @importFrom PRIMME svds
#' @importFrom irlba irlba
#' @export
eclairs_sq = function(Sigma.eclairs, rank1 = Sigma.eclairs$k, rank2=Inf, varianceFraction=1, warmStart=TRUE){

	if( ! is(Sigma.eclairs, 'eclairs') ){
		stop("ecl must be of class eclairs")
	}

	d = sqrt(Sigma.eclairs$dSq)
	V = Sigma.eclairs$U

	rank1 = min(rank1, ncol(V))

	if( rank1 < ncol(V) ){		
		V = V[,seq_len(rank1),drop=FALSE]
		d = d[seq_len(rank1),drop=FALSE]
	}

	# compute scaled eigen vectors
	# mat = sweep(V, 2, d, FUN='*')
	mat = eachrow(V, d, '*')
	rm(V)

	# compute matrix of element-wise products for all pairs of columns
	n <- ncol(mat)
	j1 <- rep.int(seq(1,n), seq(n,1))
	j2 <- sequence(seq(n,1)) - 1L + j1
	G = mat[, j1,drop=FALSE] * mat[, j2,drop=FALSE]
	rm(mat)

	# scale if j1 != j2		
	s = (j1 == j2)
	s[s] = 1
	s[!s] = sqrt(2)
	# G = sweep(G, 2, c(s)*sqrt(2), FUN='*')	
	# G = eachrow(G, c(s)*sqrt(2), '*')	
	G = eachrow(G, c(s), '*')

	# retain columns with highest variance
	####################################

	# due to precision isues, make sure variance >= 0
	df = data.frame(index = seq_len(ncol(G)), 
					variance = pmax(colVars(G), 0))

	df = df[order(df$variance, decreasing=TRUE),]
	df$cumsum = cumsum(df$variance/sum(df$variance)) 

	# set cutoff to the specified varianceFraction
	# if no columns pass this cutoff, use new cutoff
	varianceFraction = max(varianceFraction, df$cumsum[min(nrow(df), 10)])

	idx = which(df$cumsum <= varianceFraction )
	G = G[,df$index[idx],drop=FALSE]
	
	# SVD of G
	############

	n = Sigma.eclairs$n
	p = Sigma.eclairs$p
	nu = 1

	# SVD of X to get low rank estimate of Sigma
	if( rank2 < min(p, n)/3){
		if( is.null(warmStart) ){
			# dcmp = svds(G, k, isreal=TRUE)
			dcmp = irlba(G, rank2) # should be faster thatn PRIMME::svds
		}else{			
			dcmp = svds(G, rank2, u0=Sigma.eclairs$U[,seq_len(rank2)], isreal=TRUE)
		}
	}else{
		dcmp = svd(G) 

		# if rank2 < min(n,p) truncate spectrum
		if( rank2 < length(dcmp$d)){
			dcmp$u = dcmp$u[,seq_len(rank2), drop=FALSE]
			dcmp$v = dcmp$v[,seq_len(rank2), drop=FALSE]
			dcmp$d = dcmp$d[seq_len(rank2)]
		}
	}

	# eclairs
	#########

	# Estimate lambda by empirical Bayes, using nu as scale of target
	# Since data is scaled to have var 1 (instead of n), multiply by n
	res = estimate_lambda_eb( n*dcmp$d^2, n, p, nu)

	# Modify sign of dcmp$v and dcmp$u so principal components are consistant
	# This is motivated by whitening:::makePositivDiagonal()
	# but here adjust both U and V so reconstructed data is correct
	values = sign(diag(dcmp$u))

	# faster version
	dcmp$v = eachrow(dcmp$v, values, "*")
	dcmp$u = eachrow(dcmp$u, values, "*")

	ecl = list(	U 		= dcmp$u, 
				dSq 	= dcmp$d^2, 
				V 		= dcmp$v,
				lambda 	= res$lambda,
				logML 	= res$logML,
				nu 		= nu,
				n 		= n,
				p 		= p,
				k 		= length(dcmp$d),
				rownames= Sigma.eclairs$rownames,
				colnames= Sigma.eclairs$colnames,
				method 	= "svd",
				call 	= match.call())

	new("eclairs",	ecl)
}





# #' eclairs decomposition for cov/cor of transformed data
# #' 
# #' Estimate eclairs for cor/cov between transformed variables, given the SVD of the correlation matrix between original variables.  Uses parametric bootstrap
# #'
# #' @param decomp spectralDecomp 
# #' @param mu mean of each feature
# #' @param fxn function specifying transformation of original variables
# #' @param n_boot number of parametric bootstrap samples
# #' @param maxRank rank of SVD used to summarize correlation
# #' @param seed seed for parametric bootstrap sampling
# #'
# #' @return spectralDecomp object
# #'
# #' @importFrom Rfast matrnorm
# #' @importFrom RcppZiggurat zsetseed
# #' @importFrom irlba irlba
# #'
# eclairs_fxn = function( ecl, mu, fxn = identity, n_boot = 1000, maxRank=Inf, seed=1 ){

# 	if( ! is(decomp, 'spectralDecomp') ){
# 		stop("decomp must be of class spectralDecomp")
# 	}
# 	if( missing(mu) ){
# 		mu = 0
# 	}

# 	# Use SVD of covariance matrix to create multivariate normal samples
# 	A = sweep(decomp@vectors, 2, decomp@evalues^.5, FUN='*')

# 	p = nrow(decomp@vectors)

# 	# set seed for matrnorm
# 	zsetseed(seed)

# 	memUsage = p * n_boot * 8 / 2^20

# 	if( memUsage > 3000){
# 		message("Parametric bootstrap using at least ", 	
# 			format(memUsage, digits=0), " Mb memory")
# 	}

# 	# Draw parametric bootstrap samples from this covariance matrix
# 	Y_boot = decomp@vectors %*% crossprod(A, matrnorm(p, n_boot)) + mu

# 	# Compute covariance between features using transform
# 	# return SVD of the correlation matrix
# 	n = ncol(Y_boot)
# 	Y_std = scale(t(fxn(Y_boot))) / sqrt(n - 1)
# 	# Y_std = t(fxn(Y_boot))

# 	if( maxRank > min(p, n_boot)/2 ){
# 		dcmp2 <- svd( Y_std, nu=0 )
# 	}else{
# 		dcmp2 <- irlba( Y_std, nu=0, nv=maxRank)
# 	}

# 	dcmpReturn = new("spectralDecomp", 	
# 							vectors 		= dcmp2$v, 
# 							evalues 		= dcmp2$d^2, 
# 						    n_samples 		= ncol(Y_boot), 
# 						    n_features 		= nrow(Y_boot), 
# 						    varianceRetained = sum(dcmp2$d^2) / p)

# 	setDecompRank( dcmpReturn, maxRank)
# }











