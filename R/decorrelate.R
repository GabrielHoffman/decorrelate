# Gabriel Hoffman
# March 4, 2021
#
# Linear time decorrelation projection using eclairs estimate of correlation 




#' @rdname decorrelate
#' @export
setGeneric("decorrelate_eclairs", function(x, cor.est) standardGeneric("decorrelate_eclairs"))


#' @rdname decorrelate
#' @export
setMethod('decorrelate_eclairs', c(	x = "matrix", cor.est = "eclairs"), 
	function(x, cor.est) {

	decorrelate(x, cor.est$U, cor.est$dSq*(1- cor.est$lambda), cor.est$lambda)
})



#' Decorrelation projection
#' 
#' Decorrelation projection using low rank singular value decomposition of the correlation matrix
#'
#' @param x matrix to be transformed so *columns* are independent
#' @param U orthonormal matrix with k columns representing the low rank component
#' @param dSq eigen values so that U diag(lamba) U^T is the low rank component
#' @param sigSq residual variance used for the scaled diagonal component: diag(sigSq, p)
#' @param cor.est estimate of correlation matrix from \link{eclairs} storing \code{U}, \code{dSq}, and \code{sigSq}
#'
#' @description
#' Given a vector \eqn{x ~ N(0, \Sigma)} where \eqn{\Sigma = U diag(\lambda) U^T + diag(\sigma^2)}, transform x so that it has an identity covariance matrix.  \eqn{\Sigma^{-0.5}x} is such a projection (see Strimmer whitening).  When \eqn{\Sigma} is \eqn{p \times p}, computing this project naively is \eqn{O(p^3)}.  Here we take advantage of the fact that \eqn{\Sigma} is the sum of a low rank decomposition, plus a diagonal matrix with a single value.
#'
#' @rdname decorrelate
#' @export
setGeneric("decorrelate", function(x, U, dSq, sigSq) standardGeneric("decorrelate"))



#' @rdname decorrelate
#' @export
setMethod('decorrelate', signature = c(
				x 		= "matrix", 
				U 		= "matrix", 
				dSq 	= "numeric", 
				sigSq	= "numeric"), 
	function(x, U, dSq, sigSq) {

	transpose = FALSE
	if( min(dim(x)) == 1 ){
		transpose = TRUE
		x = matrix(x, nrow=1)
	}

	# # add dimension checks here
	if( ncol(x) != nrow(U) ){
		stop("length of x must equal nrow(U)")
	}

	if( ncol(U) != length(dSq) ){
		stop("ncol(U) must equal length of dSq")
	}

	if( length(sigSq) != 1){
		stop( "length of sigSq must be 1")
	}

	x.decorr = lrmult(x, U, dSq, sigSq, alpha = -1/2)

	if( transpose ) x.decorr = t(x.decorr)

	x.decorr
})

# low rank multiplication
# x %*% (U diag(dSq) U^T + diag(sigSq))^alpha
lrmult = function(x, U, dSq, sigSq, alpha){

	# precompute value that is used twice
	a = (x - tcrossprod(x %*% U, U)) * sigSq^alpha

	tcrossprod(x %*% U * (dSq + sigSq)^alpha, U) + (a - tcrossprod(a %*% U, U))
}


lrmult_t = function(x, U, dSq, sigSq, alpha){
	lrmult(t(x), U, dSq, sigSq, alpha)
}



#  # low rank multiplication
# # (U diag(dSq) U^T + diag(sigSq))^alpha %*% x
lrmult_right = function(x, U, dSq, sigSq, alpha){

	# precompute value that is used twice
	a = (x - U %*% crossprod(U, x)) * sigSq^alpha

	U %*% (crossprod(U, x) * (dSq + sigSq)^alpha) + (a - U %*% crossprod(U, a))
}

lrmult_left = function(x1, U, dSq, sigSq, alpha){

	# precompute value that is used twice
	b = (x1 - tcrossprod(x1 %*% U, U)) * sigSq^alpha

	# (x1 %*% U * (dSq + sigSq)^alpha) %*% U + (b - tcrossprod(b,U) %*% U)
	(x1 %*% U) %*% (( (dSq + sigSq)^alpha) * t(U)) + (b - tcrossprod(b %*% U, U))
}





# x = t(Y)
# x1 = Y

# r1 = (a - U %*% crossprod(U, a))
# r2 = (b - tcrossprod(b %*% U, U))
# range(r1 - t(r2))





# X = scale(Y) / sqrt(n-1)
# res = ShrinkCovMat::shrinkcovmat.unequal( t(X) )
# str(res)

# Sigma = tcrossprod(X) 

# dcmp = eigen(Sigma)
# ev = dcmp$values

# tr = function(C) sum(diag(C))

# tr(Sigma)



# (n+n^2) / (n^2 + (p-n+1)/p * n^2 )


# ShrinkCovMat:::trace_stats_centered((X))




# stimate lambda two ways
# res = ShrinkCovMat::shrinkcovmat.equal(t(X)* sqrt(n-1))
# abs(res$lambdahat - est.shrink$lambda)

# est.shrink = mvIC:::shrinkcovmat.equal_lambda(t(X)* sqrt(n-1))
# Sigma = with(dcmp, (v %*% diag(d^2) %*% t(v)))
# Sig.hat = with(est.shrink, (1-lambda_hat) * Sigma + diag(lambda_hat*1, p))

# range(res$Sigmasample - Sigma)
# range(res$Sigmahat - Sig.hat)

# Sig.hat[1:3, 1:3]
# res$Sigmahat[1:3, 1:3]



test_decorrelate = function(){

	set.seed(1)
	n = 20 # number of samples
	p = 500 # number of features
	k = 2 # rank of covariance matrix
	lambda = 0.1

	U_latent = svd(matrnorm(p,100))$u[,1:k, drop=FALSE]
	Sigma = diag(.1, p) + tcrossprod(U_latent)

	# sample matrix from MVN with covariance Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)

	# sample vector from MVN with covariance Sigma
	zstat = t(rmvnorm(1, rep(0, p), sigma=Sigma))

	# Estimate covariance from Y
	cor.est = eclairs( Y, k=3, lambda=lambda)

	# perform low rank decorrelation projection
	z.decorr = decorrelate((zstat), cor.est$U, cor.est$dSq*(1- cor.est$lambda), cor.est$lambda)

	# perform *full* rank decorrelation projection
	edcomp = eigen(Sigma)

	# Decorrelation using full spectrum from eigen decomp
	M.transform = with(edcomp, vectors %*% diag(1/sqrt(values + lambda)) %*% t(vectors) + diag(1/sqrt(lambda), p))
	z.transform = M.transform %*% zstat

	cor(z.transform, z.decorr)



	cor.est = eclairs( Y, k=2, lambda)


	z.decorr = decorrelate(zstat, cor.est$U, cor.est$dSq *(1-cor.est$lambda), cor.est$lambda)

	cor(z.transform, z.decorr)



	# NO estimation here, just the linear algebra part
	#######
	dcmp = irlba::partial_eigen(Sigma, 20)
	uhat = diag(Sigma - with(dcmp, vectors %*% diag(values) %*% t(vectors)))

	# Compute residual variance of each feature
	sigSq = mean(uhat)

	z.decorr = decorrelate(zstat, dcmp$vectors, dcmp$values, sigSq)

	cor(z.transform, z.decorr)
	
}













