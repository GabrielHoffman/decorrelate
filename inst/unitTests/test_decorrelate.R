
library(RUnit)

test_decorrelate = function(){

	suppressPackageStartupMessages({
	library(Rfast)
	})

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


	z.decorr = decorrelate(zstat, U=cor.est$U, dSq=cor.est$dSq*(1- cor.est$lambda), sigSq=cor.est$lambda)


	z.decorr2 = decorrelate_eclairs(zstat, cor.est )

	checkIdentical(z.decorr, z.decorr2)
}


test_lrmult = function(){

	suppressPackageStartupMessages({
	library(Rfast)
	})

	set.seed(1)
	n = 20 # number of samples
	p = 500 # number of features
	k = 2 # rank of covariance matrix
	lambda = 0.1

	U_latent = svd(matrnorm(p,100))$u[,1:k, drop=FALSE]
	Sigma = diag(.1, p) + tcrossprod(U_latent)

	# sample matrix from MVN with covariance Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)



	dcmp = eigen(cov2cor(Sigma))
	k = 10
	U = dcmp$vectors
	U1 = U[,1:k]
	U2 = U[,-c(1:k)]
	dcmp$values[-c(1:k)] = 0
	dSq = dcmp$values
	dSq1 = dSq[1:k]
	lambda = 0.1

	C1 = U %*% diag(dSq) %*% t(U) * (1-lambda) + diag(lambda, p)
	# C2 = U %*% diag(dSq*(1-lambda)) %*% t(U) + diag(lambda, p)
	# C2 = U %*% diag(dSq*(1-lambda) + lambda) %*% t(U) 
	# C2 = U %*% diag(dSq*(1-lambda) + lambda) %*% t(U) 

	v = dSq[1:k]*(1-lambda) + lambda
	# C2 = U1 %*% diag(v) %*% t(U1) + U2 %*% diag(lambda, p-k) %*% t(U2)
	# C2 = U1 %*% diag(v) %*% t(U1) + tcrossprod(U2) *lambda
	# C2 = U1 %*% (v * t(U1)) + tcrossprod(U2) *lambda
	C2 = U1 %*% (v * t(U1)) + tcrossprod(diag(1,p) - tcrossprod(U1)) *lambda

	checkEquals(C1, C2)


	# projection
	alpha = -1/2
	C = U %*% diag(dSq) %*% t(U) * (1-lambda) + diag(lambda, p)

	decomp = eigen(C)
	Y.decorr1 = Y %*% with(decomp, vectors %*% diag(values^alpha) %*% t(vectors))

	# Y.decorr2 = Y %*% (U1 %*% (v^alpha * t(U1)) + tcrossprod(diag(1,p) - tcrossprod(U1)) *lambda^alpha)

	# Y.decorr2 = (Y %*% U1) %*% (v^alpha * t(U1)) + Y %*%tcrossprod(diag(1,p) -  tcrossprod(U1)) *(lambda^alpha)

	# Y.decorr2 = (Y %*% U1) %*% (v^alpha * t(U1)) +Y %*% (diag(1,p)  - tcrossprod(U1) ) *(lambda^alpha)

	Y.decorr2 = (Y %*% U1) %*% (v^alpha * t(U1)) + (Y  - tcrossprod(Y %*% U1,U1) ) *(lambda^alpha)

	checkEquals(Y.decorr1, Y.decorr2)


	Y.decorr2 = mult_eclairs(Y, U1, dSq1, lambda, alpha)

	checkEquals(Y.decorr1, Y.decorr2)
}


test_decorrelate_full_rank = function(){
	# Check that decorrelate gives same value as 1/sqrt eigen values
	######################################
	set.seed(1)
	n = 400 # number of samples
	p = 500 # number of features
	k = 2 # rank of covariance matrix
	lambda = 0.1

	# Create correlation matrix with autocorrelation
	autocorr.mat <- function(p = 100, rho = 0.9) {
	    mat <- diag(p)
	    return(rho^abs(row(mat)-col(mat)))
	}

	Sigma = autocorr.mat(p, .9)
	Sigma[1:3, 1:3]

	# sample matrix from MVN with covariance Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)
	# cor(Y[,1:3])


	# decorrelate based on true correlation matrix
	Y.decorr = Y %*% with(eigen(Sigma), vectors %*% diag(values^-0.5) %*% t(vectors))
	# cor(Y.decorr[,1:3])

	dcmp = eigen(Sigma)
	cor.est = list(U = dcmp$vectors,
				dSq = dcmp$values,
				lambda =1e-10)

	Y.decorr2 = decorrelate(Y, cor.est)

	# Y.decorr2 = mult_eclairs(Y, cor.est$U, cor.est$dSq, cor.est$lambda, -1/2)
	cor(Y.decorr2[,1:3])

	checkEquals(Y.decorr, Y.decorr2)

}

test = function(){

	# test statistical performance of decorrelate 
	# 	It performs well with n >> p, but poorly in p >> n
	#   However, this is not an issue because in regression, 
	# 	the covariance matrix is actually known

	library(Matrix)
	set.seed(1)
	n = 100 # number of samples
	p = 2000 # number of features

	# Create correlation matrix with autocorrelation
	autocorr.mat <- function(p = 100, rho = 0.9) {
	    mat <- diag(p)
	    return(rho^abs(row(mat)-col(mat)))
	}

	Sigma = autocorr.mat(p/8, .9999)
	Sigma = bdiag(Sigma, Sigma)
	Sigma = bdiag(Sigma, Sigma)
	Sigma = bdiag(Sigma, Sigma)

	Sigma[1:3, 1:3]

	# sample matrix from MVN with covariance Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)
	# cor(Y[,1:3])

	# decorrelated using true correlation matrix
	Y.decorr.exact = Y %*% with(eigen(Sigma), vectors %*% diag(values^-0.5) %*% t(vectors))

	cor(Y.decorr.exact[,1:3])
	hist(cor(Y.decorr.exact))

	ecl = eclairs(Y)

	plot(ecl$dSq)

	C = getCor(ecl)

	Sigma[1:3, 1:3]
	C[1:3, 1:3]


	Y.decorr2 = decorrelate(Y, ecl)
	cor(Y[,1:3])
	cor(Y.decorr2[,1:3])



	hist(cor(Y))
	hist(cor(Y.decorr2))
}























