
test_new = function(){

	library(decorrelate)
	library(Rfast)
	n = 8000
	p = 2

	Sigma = autocorr.mat(p, .5)*2

	Y = rmvnorm(n, rep(0, p), sigma=Sigma, seed=1)

	# whitening transform
	A = decorrelate::whiten(Y, lambda=0)
	B = whitening::whiten(Y)
	checkEqualsNumeric(cov(A), cov(B))

	# applying whitening matrix directly
	ecl = eclairs(Y, lambda=0, compute="cov")
	W1 = getWhiteningMatrix(ecl, lambda=0)
	A = tcrossprod(Y, W1)
	W2 = whitening::whiteningMatrix(cov(Y))
	B = tcrossprod(Y, W2)
	checkEqualsNumeric(cov(A), cov(B))
}

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

	z.decorr = t(mult_eclairs(t(zstat), cor.est$U, cor.est$dSq, cor.est$lambda, cor.est$nu, alpha = -1/2, cor.est$sigma))

	z.decorr2 = decorrelate(zstat, cor.est, transpose=TRUE )

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

	a = checkEquals(C1, C2)


	# projection
	alpha = -1/2
	C = U %*% diag(dSq) %*% t(U) * (1-lambda) + diag(lambda, p)

	decomp = eigen(C)
	Y.decorr1 = Y %*% with(decomp, vectors %*% diag(values^alpha) %*% t(vectors))

	# Y.decorr2 = Y %*% (U1 %*% (v^alpha * t(U1)) + tcrossprod(diag(1,p) - tcrossprod(U1)) *lambda^alpha)

	# Y.decorr2 = (Y %*% U1) %*% (v^alpha * t(U1)) + Y %*%tcrossprod(diag(1,p) -  tcrossprod(U1)) *(lambda^alpha)

	# Y.decorr2 = (Y %*% U1) %*% (v^alpha * t(U1)) +Y %*% (diag(1,p)  - tcrossprod(U1) ) *(lambda^alpha)

	Y.decorr2 = (Y %*% U1) %*% (v^alpha * t(U1)) + (Y  - tcrossprod(Y %*% U1,U1) ) *(lambda^alpha)

	b = checkEquals(Y.decorr1, Y.decorr2)

	sig = rep(1, ncol(Y))	
	Y.decorr2 = mult_eclairs(Y, U1, dSq1, lambda, 1, alpha, sig)

	c = checkEquals(Y.decorr1, Y.decorr2)

	a & b & c
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
				lambda =1e-10,
				nu=1)
	cor.est = new("eclairs", cor.est)

	Y.decorr2 = decorrelate(Y, cor.est)

	# Y.decorr2 = mult_eclairs(Y, cor.est$U, cor.est$dSq, cor.est$lambda, -1/2)
	cor(Y.decorr2[,1:3])

	checkEquals(Y.decorr, Y.decorr2)
}

# test = function(){

# 	# test statistical performance of decorrelate 
# 	# 	It performs well with n >> p, but poorly in p >> n
# 	#   However, this is not an issue because in regression, 
# 	# 	the covariance matrix is actually known

# 	library(Matrix)
# 	set.seed(1)
# 	n = 100 # number of samples
# 	p = 2000 # number of features

# 	# Create correlation matrix with autocorrelation
# 	autocorr.mat <- function(p = 100, rho = 0.9) {
# 	    mat <- diag(p)
# 	    return(rho^abs(row(mat)-col(mat)))
# 	}

# 	Sigma = autocorr.mat(p/8, .9999)
# 	Sigma = bdiag(Sigma, Sigma)
# 	Sigma = bdiag(Sigma, Sigma)
# 	Sigma = bdiag(Sigma, Sigma)

# 	Sigma[1:3, 1:3]

# 	# sample matrix from MVN with covariance Sigma
# 	Y = rmvnorm(n, rep(0, p), sigma=Sigma)
# 	# cor(Y[,1:3])

# 	# decorrelated using true correlation matrix
# 	Y.decorr.exact = Y %*% with(eigen(Sigma), vectors %*% diag(values^-0.5) %*% t(vectors))

# 	cor(Y.decorr.exact[,1:3])
# 	hist(cor(Y.decorr.exact))

# 	ecl = eclairs(Y)

# 	plot(ecl$dSq)

# 	C = getCor(ecl)

# 	Sigma[1:3, 1:3]
# 	C[1:3, 1:3]


# 	Y.decorr2 = decorrelate(Y, ecl)
# 	cor(Y[,1:3])
# 	cor(Y.decorr2[,1:3])



# 	hist(cor(Y))
# 	hist(cor(Y.decorr2))
# }










test_lm_each_eclairs = function(){

	library(Matrix)
	library(Rfast)
	set.seed(1)
	n = 800 # number of samples
	p = 8*200 # number of features

	# Create correlation matrix with autocorrelation
	autocorr.mat <- function(p = 100, rho = 0.9) {
	mat <- diag(p)
	return(rho^abs(row(mat)-col(mat)))
	}

	# create correlation matrix
	Sigma = autocorr.mat(p/8, .9)
	Sigma = bdiag(Sigma, Sigma)
	Sigma = bdiag(Sigma, Sigma)
	Sigma = bdiag(Sigma, Sigma)

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma*5.1)

	# eclairs decomposition
	ecl = eclairs(Y)

	# simulate covariates
	data = data.frame(matrnorm(p,2))
	colnames(data) = paste0('v', 1:2)

	# simulate response
	y = rnorm(p)

	# Simulate 1000 features to test
	X = matrnorm(p, 1000)
	colnames(X) = paste0('set_', 1:ncol(X))

	lm_each_eclairs_slow = function(formula, data, X, Sigma.eclairs,...){

	    formula.mod = update(formula, . ~ . + x.tempVariable)

	    res = t(apply(X, 2, function(x){

	        # included this x vector in the data matrix
	        data$x.tempVariable = x

	        # fit regression
	        fit = lm_eclairs( formula.mod, data, Sigma.eclairs,...)

	        # extract hypothesis test
	        coef(summary(fit))['x.tempVariable',,drop=FALSE]    
	    }))

	    rownames(res) = colnames(X)
	    colnames(res) = c("beta", "se", "tstat", "pvalue")

	    as.data.frame(res)
	}

	# Test time

	# # Use slow version
	# system.time({
	# res1 <- lm_each_eclairs_slow(y ~ v1 + v2, data, X[,1:1000], ecl )
	# })

	# # fast version using pre-project ( ~ 50X faster!)
	# system.time({
	# res2 <- lm_each_eclairs(y ~ v1 + v2, data, X[,1:1000], ecl )
	# })

	# Use slow version
	res1 = lm_each_eclairs_slow(y ~ v1 + v2, data, X[,1:10], ecl )

	# fast version using pre-project
	res2 = lm_each_eclairs(y ~ v1 + v2, data, X[,1:10], ecl )

	checkEquals(res1, res2)
}


test_reform_decomp = function(){

	library(RUnit)
	# library(PRIMME)
	library(decorrelate)
	library(Matrix)
	library(Rfast)
	set.seed(1)
	n = 800 # number of samples
	p = 8*200 # number of features

	# Create correlation matrix with autocorrelation
	autocorr.mat <- function(p = 100, rho = 0.9) {
	mat <- diag(p)
	return(rho^abs(row(mat)-col(mat)))
	}

	# create correlation matrix
	Sigma = autocorr.mat(p/8, .9)
	Sigma = bdiag(Sigma, Sigma)
	Sigma = bdiag(Sigma, Sigma)
	Sigma = bdiag(Sigma, Sigma)

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma*5.1)
	rownames(Y) = paste0("sample_", 1:n)
	colnames(Y) = paste0("gene_", 1:p)

	# Correlation
	#------------

	# eclairs decomposition
	Sigma.eclairs = eclairs(Y, compute="correlation", lambda=.1)

	# Compute SVD on subset of eclairs decomposition
	drop = paste0("gene_",1:200)

	# reform
	ecl1 = reform_decomp( Sigma.eclairs, drop=drop)

	# eclairs using original data
	ecl2 = eclairs(Y[,!colnames(Y) %in% drop], compute="correlation", lambda=.1)
	

	# since the last singular value is very small (below machine precision)
	# the large singular vectors of U and V are unstable 
	# and should be ommitted from the check
	idx = seq(1:(Sigma.eclairs$k-1))
	# max(ecl1$U[,idx] - ecl2$U[,idx])
	# max(ecl1$V[,idx] - ecl2$V[,idx])

	res1 = checkEqualsNumeric( ecl1$V[,idx], ecl2$V[,idx])
	res2 = checkEqualsNumeric( ecl1$U[,idx], ecl2$U[,idx])

	# Covariance
	############# 

	# eclairs decomposition
	Sigma.eclairs = eclairs(Y, compute="covariance", lambda=.1)

	# Compute SVD on subset of eclairs decomposition
	drop = paste0("gene_",1:200)

	# reform
	ecl1 = reform_decomp( Sigma.eclairs, drop=drop)

	# eclairs using original data
	ecl2 = eclairs(Y[,!colnames(Y) %in% drop], compute="covariance", lambda=.1)
	

	# since the last singular value is very small (below machine precision)
	# the large singular vectors of U and V are unstable 
	# and should be ommitted from the check
	idx = seq(1:(Sigma.eclairs$k-1))
	# max(ecl1$U[,idx] - ecl2$U[,idx])
	# max(ecl1$V[,idx] - ecl2$V[,idx])

	res3 = checkEqualsNumeric( ecl1$V[,idx], ecl2$V[,idx])
	res4 = checkEqualsNumeric( ecl1$U[,idx], ecl2$U[,idx])

	res1 & res2 & res3 & res4
}





test_decorrelate_transpose = function(){

	library(Matrix)
	library(Rfast)
	set.seed(1)
	n = 800 # number of samples
	p = 8*200 # number of features

	# Create correlation matrix with autocorrelation
	autocorr.mat <- function(p = 100, rho = 0.9) {
	mat <- diag(p)
	return(rho^abs(row(mat)-col(mat)))
	}

	# create correlation matrix
	Sigma = autocorr.mat(p/8, .9)
	Sigma = bdiag(Sigma, Sigma)
	Sigma = bdiag(Sigma, Sigma)
	Sigma = bdiag(Sigma, Sigma)

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma*5.1)
	

	# eclairs decomposition
	ecl = eclairs(Y, k=10)

	a <- decorrelate(Y, ecl)
	b <- t(decorrelate(t(Y), ecl, transpose=TRUE))

	checkEqualsNumeric(a,b)
}



test_whitening_matrix = function(){

	library(Matrix)
	library(Rfast)
	library(ggplot2)
	library(cowplot)
	library(whitening)

	n = 2000
	p = 3

	Y = matrnorm(n,p)*10

	# decorrelate with implicit whitening matrix
	# give same result as explicity whitening matrix
	ecl <- eclairs(Y, compute="covariance", lambda=0)

	W = getWhiteningMatrix( ecl, lambda=0 )
	Z1 = tcrossprod(Y, W)
	Z2 = decorrelate(Y, ecl)
	checkEqualsNumeric(Z1, Z2)

	# compare decorrelate and whiten
	ecl <- eclairs(Y, compute="covariance", lambda=0)
	Z1 <- decorrelate(Y, ecl, lambda=0)
	Z2 = decorrelate::whiten(Y, lambda=0)
	checkEqualsNumeric(Z1, Z2)

	# compare covariance matrices of whitened data
	ecl <- eclairs(Y, compute="covariance", lambda=0)
	Z1 <- decorrelate(Y, ecl, lambda=0)
	Z2 = whitening::whiten(Y)
	checkEqualsNumeric(cov(Z1), cov(Z2))

	# covariance
	###############
	ecl <- eclairs(Y, compute="covariance")

	# compare whitening matrices
	W1 = getWhiteningMatrix( ecl, lambda=0 )
	W2 = whiteningMatrix(cov(Y), method='ZCA')
	checkEqualsNumeric(cov(tcrossprod(Y, W1)), cov(tcrossprod(Y, W2)))

	# compare transformed data
	Z1 = tcrossprod(Y, W1)
	Z2 = tcrossprod(Y, W2)
	checkEqualsNumeric(cov(Z1), cov(Z2))

	# since covariance is removed, 
	# the covariance of the trasnformed data is identity
	# cov(Z1)
	# cor(Z1)

	Z3 = whitening::whiten(Y, method="ZCA")
	checkEqualsNumeric(Z3, Z2)


	# correlation
	ecl <- eclairs(Y, compute="correlation")

	# compare whitening matrices
	W1 = getWhiteningMatrix( ecl, lambda=0 )
	W2 = whiteningMatrix(cor(Y), method='ZCA')
	checkEqualsNumeric(W1, W2)

	# compare transformed data
	Z1 = tcrossprod(Y, W1)
	Z2 = tcrossprod(Y, W2)
	checkEqualsNumeric(Z1, Z2)


	# since correlation is removed, 
	# the correlation is shrunk toward zero.
	# But because the diagonal cov is not exactly zero, the 
	# decorrelation is approximate.  
	# But if diagonal elements are very similar, method is very good
	# are 
	# cov(Z1)
	# cor(Z1)

	# Z3 = whiten(scale(Y), method="PCA")
	# checkEqualsNumeric(Z3, Z2)

}



# library(RUnit)


# set.seed(1)
# n = 300000 # number of samples
# p = 2 # number of features]

# # sample matrix from MVN with covariance Sigma
# Sigma = matrix(c(3, 2, 2, 5), 2,2)
# Y = rmvnorm(n, rep(0, p), sigma=Sigma, set.seed=1)

# # same
# ecl = eclairs(Y, compute = "cor", lambda=0)
# ecl$U
# ecl$dSq
# eigen(cor(Y))


# # why is this different?
# # because sigma is variable
# ecl = eclairs(Y, compute = "cov", lambda=0)
# ecl$U
# ecl$dSq
# eigen(cov(Y))





# minvsqrt <- function(S) {
# dcmp <- eigen(S)
# with(dcmp, vectors %*% diag(1 / sqrt(values)) %*% t(vectors))
# }
# cov(Y %*% minvsqrt(cor(Y)))

# ecl = eclairs(Y, compute = "cor", lambda=0)
# getCor(ecl)
# cor(decorrelate(Y, ecl))
# cor(Y)
# cor(Y %*% minvsqrt(cov(Y)))



# ecl = eclairs(Y, compute = "cor", lambda=0)
# whitening::whiteningMatrix(cor(Y))
# getWhiteningMatrix(ecl)

# # decorrelate doesn't set cor to zero
# cov(whitening::whiten(Y))
# cov(decorrelate::whiten(Y, lambda=0))



# ecl = eclairs(Y, compute = "cov", lambda=0)
# whitening::whiteningMatrix(cov(Y))
# getWhiteningMatrix(ecl)


# getWhiteningMatrix(ecl, lambda=0)
# whitening::whiteningMatrix(cov(Y))




# W = getWhiteningMatrix(ecl)
# cor(tcrossprod(Y, W))
# cor(whiten(Y, lambda=0))


# W = whitening::whiteningMatrix(cov(Y))
# cor(tcrossprod(Y, W))



# ecl = eclairs(Y, compute = "cov", lambda=0)
# getCov(ecl)
# cov(Y)
# cov(decorrelate(Y, ecl))



# W = getWhiteningMatrix(ecl)
# whitening::whiteningMatrix(cov(Y))
# cov(tcrossprod(Y, W))

# # whitening matrices are the same
# ecl = eclairs(Y, compute = "cor", lambda=0)
# getWhiteningMatrix(ecl)
# minvsqrt(cor(Y))

# # slight rounding error?
# # or is factoring out Z wrong?
# ecl = eclairs(Y, compute = "cov", lambda=0)
# getWhiteningMatrix(ecl)
# minvsqrt(cov(Y))


# cov(Y %*% getWhiteningMatrix(ecl))
# cov(Y %*% minvsqrt(cov(Y)))


# ecl = eclairs(Y, compute = "cov", lambda=0)
# v <- ecl$dSq

# alpha <- -1 / 2

# if( any(ecl$sigma != 1)){
# # ecl$U <- diag(ecl$sigma^alpha) %*% ecl$U 
# ecl$U <- diag(1/sqrt(ecl$sigma)) %*% ecl$U 
# }

# W <- with(ecl, U %*% diag(1/sqrt(v)) %*% t(U))



















