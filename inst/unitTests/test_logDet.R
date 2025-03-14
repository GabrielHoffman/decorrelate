

test_logDet = function(){

	library(Rfast)
	set.seed(1)
	n = 800 
	p = 12

	# create correlation matrix
	Sigma = autocorr.mat(p, .9) *4

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)

	# eclairs decomposition
	ecl = eclairs( Y )

	a = logDet(ecl)
	b = determinant(getCov(ecl))$modulus
	checkEqualsNumeric(a,b)

	a = logDet(ecl, 0.5)
	b = determinant(chol(getCov(ecl)))$modulus
	checkEqualsNumeric(a,b)

	a = logDet(ecl, -1)
	b = determinant(solve(getCov(ecl)))$modulus
	checkEqualsNumeric(a,b)
}


test_tr = function(){

	library(Rfast)
	set.seed(1)
	n = 800 
	p = 12

	# create correlation matrix
	Sigma = autocorr.mat(p, .9)

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)

	# eclairs decomposition
	ecl = eclairs( Y, compute="corr" )

	a = tr(ecl)
	b = sum(eigen(getCor(ecl))$values)
	d = sum(diag(getCor(ecl)))
	checkEqualsNumeric(a,b)
	checkEqualsNumeric(a,d)

	# test for covariance with non-identity sigma
	# create covariance matrix
	# Sigma = autocorr.mat(p, .9) * 4
	# Y = rmvnorm(n, rep(0, p), sigma=Sigma)

	# # eclairs decomposition
	# ecl = eclairs( Y )

	# a = tr(ecl)
	# b = sum(eigen(getCov(ecl))$values)
	# d = sum(diag(getCov(ecl)))
	# checkEqualsNumeric(a,b)
	# checkEqualsNumeric(a,d)

	# a = diag(getCor(ecl))
	# sum(a*ecl$sigma^2)

	# a = eigen(getCov(ecl))$values
	# sum(a*ecl$sigma^2)


}



test_kappa = function(){

	library(Rfast)
	set.seed(1)
	n = 800 
	p = 12

	# create correlation matrix
	Sigma = autocorr.mat(p, .9)

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)

	# eclairs decomposition
	Sigma.eclairs = eclairs( Y, compute="corr" )

	a = kappa(Sigma.eclairs)
	b = kappa(getCor(Sigma.eclairs), exact=TRUE)
	checkEqualsNumeric(a,b)
}









test_quadForm = function(){

	library(Rfast)
	set.seed(1)
	n = 800 
	p = 12

	# create covariance matrix
	Sigma = autocorr.mat(p, .9) * 4

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma*4.5)

	# eclairs decomposition
	ecl = eclairs( Y )

	# return scalar
	a = quadForm( ecl, Y[1,])
	b = Y[1,] %*% solve(getCov(ecl)) %*% ( Y[1,])
	checkEqualsNumeric(a,b)

	# return matrix
	a = quadForm( ecl, Y[1:2,])
	b = Y[1:2,] %*% solve(getCov(ecl)) %*% t( Y[1:2,])
	checkEqualsNumeric(a,b)
}








