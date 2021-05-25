

test_estimate_lambda = function(){

	library(decorrelate)
	# source("~/workspace/repos/decorrelate/R/eb_lambda.R")
	# devtools::reload("./")
	set.seed(1)

	n = 100
	p = 1000
	nu = 1

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

	# test stability of lambda as rank changes
    # sapply((n-20):n, function(k){
    # 	ev = ecl$dSq[1:k]
	   #  estimate_lambda_eb(n*ev, n, p, nu)
	 # })


	# x = 3
	# ev = eigs
	# ev[n:(n-x)] = mean(ev[n:(n-x)] )
	# decorrelate::estimate_lambda_eb(ev, n, p, nu=nu)


	# ev[n:(n-x)] = 0
	# decorrelate::estimate_lambda_eb(ev, n, p, nu=nu)



	# ecl = eclairs(X, compute="corr", k=90)
	# totalVar = ecl$p
	# ev = ecl$dSq
	# idx = (length(ev)+1):ecl$n
	# ev[idx] = (totalVar-sum(ev)) / length(idx)
 #   	decorrelate:: estimate_lambda_eb(n*ev, n, p, nu)


 	 # estimate_lambda_eb(n * ecl$dSq, n, ecl$p, nu)

 	 # ev = n * ecl$dSq
 	 # p =  ecl$p


 	 # estimate_lambda_eb(ev, n, ecl$p, nu)


	# ev = eigs
	# ev[p:length(ev)] = 0
	# lambda = c()
	# for(i in length(eigs):80){
	# 	ev[i] = 0
	# 	value = decorrelate::estimate_lambda_eb(ev, n, p, nu=nu)
	# 	lambda = c(lambda, value)
	# }
	# plot(lambda)

	# cumsum(eigs[1:90]/sum(eigs))

	# ev = eigs
	# ev[p:length(ev)] = 0
	# ev[86:n] = 0
	# ratio = sum(eigs)/ sum(ev)
	# decorrelate::estimate_lambda_eb(ev*ratio, n, p, nu=nu)

	# # fit@alphaOpt
	# # est
	# # decorrelate:::shrinkcovmat.equal_lambda(X)

	ecl = eclairs(X, compute="corr")
    # lambda_eclairs = estimate_lambda_eb(n*ecl$dSq, n, p, nu)

	# check two estimates of lambda
	checkEqualsNumeric( fit@alphaOpt, est, tolerance=1e-4) & checkEqualsNumeric( fit@alphaOpt, ecl$lambda, tolerance=1e-4)
}

test_large_scale = function(){

	# run eclairs on very large dataset

	q()
	R
	library(decorrelate)
	library(Matrix)
	library(Rfast)

	set.seed(1)
	n = 1200 # number of samples
	p = 100 # number of feature per block

	# Create correlation matrix with autocorrelation
	autocorr.mat <- function(p = 100, rho = 0.9) {
	 mat <- diag(p)
	 return(rho^abs(row(mat)-col(mat)))
	}

	# create correlation 
	Sigma.ch = chol(autocorr.mat(p, .9))
	
	for( i in 1:2)
		Sigma.ch = bdiag(Sigma.ch, Sigma.ch)

	ncol(Sigma.ch)

	# draw data from correlation matrix Sigma
	# Y = rmvnorm(n, rep(0, nrow(Sigma)), sigma=Sigma*5.1)
	Y = matrnorm(n, ncol(Sigma.ch)) %*% Sigma.ch
	Y = as.matrix(Y)

	ecl = eclairs(Y[,1:50000 + 100000])

	plot(ecl)
}













