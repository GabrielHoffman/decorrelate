
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

Y.decorr = decorrelate_eclairs(Y, cor.est )
