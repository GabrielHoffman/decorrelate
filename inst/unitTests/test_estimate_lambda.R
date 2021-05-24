

test_estimate_lambda = function(){

	library(decorrelate)

	set.seed(1)

	n = 100
	p = 500
	nu = 2

	autocorr.mat <- function(p = 100, rho = 0.9) {
	    mat <- diag(p)
	    return(rho^abs(row(mat)-col(mat)))
	}

	# cobj = genPositiveDefMat(p)
	# Sigma = cov2cor(cobj$Sigma)
	Sigma = autocorr.mat(p, .9)

	X = matrix(rnorm(n*p), n, p)
	# Use scaling from beam's Rcpp internals
	X = scale(X %*% chol(Sigma)) * sqrt(n / (n-1))

	# beam
	######

	fit = beam::beam(X, type="marginal", return.only="cor", D = diag(nu, p))

	# fit@alphaOpt

	# head(fit@gridAlpha)

	# version in decorrelate
	########################

	# evaluate eigen-values
	if( n > p){
		eigs = eigen(crossprod(X), symmetric=TRUE)$values
	}else{
		eigs = eigen(tcrossprod(X), symmetric=TRUE)$values
	}

	# estimate lambda value
	est = decorrelate::estimate_lambda_eb(eigs, n, p, nu=nu)


	# fit@alphaOpt
	# est
	# decorrelate:::shrinkcovmat.equal_lambda(X)

	ecl = eclairs(X, compute="corr")
    lambda_eclairs = estimate_lambda_eb(n*ecl$dSq, n, p, nu)

	# check two estimates of lambda
	checkEqualsNumeric( fit@alphaOpt, est, tolerance=1e-4) & checkEqualsNumeric( fit@alphaOpt, lambda_eclairs, tolerance=1e-4)
}




test_large_scale = function(){


	# library(decorrelate)

	# library(Matrix)
	# library(Rfast)
	# set.seed(1)
	# n = 800 # number of samples
	# p = 512*100 # number of features

	# # Create correlation matrix with autocorrelation
	# autocorr.mat <- function(p = 100, rho = 0.9) {
	#  mat <- diag(p)
	#  return(rho^abs(row(mat)-col(mat)))
	# }

	# # create correlation matrix
	# Sigma = autocorr.mat(p/512, .9)
	# Sigma = bdiag(Sigma, Sigma)
	# Sigma = bdiag(Sigma, Sigma)
	# Sigma = bdiag(Sigma, Sigma)
	# Sigma = bdiag(Sigma, Sigma)
	# Sigma = bdiag(Sigma, Sigma)
	# Sigma = bdiag(Sigma, Sigma)
	# Sigma = bdiag(Sigma, Sigma)
	# Sigma = bdiag(Sigma, Sigma)



	# # draw data from correlation matrix Sigma
	# Y = rmvnorm(n, rep(0, nrow(Sigma)), sigma=Sigma*5.1)

	# ecl = eclairs(Y)
	# plot(ecl)


}







